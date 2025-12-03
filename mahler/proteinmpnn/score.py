import numpy as np
import torch
import random
from pathlib import Path
import mdtraj as md

from .protein_mpnn_utils import tied_featurize, StructureDatasetPDB, ProteinMPNN, parse_file, parse_traj

class ScoreMPNN:

    def __init__(self,     
                 seed: int = 42,
                 ca_only: bool = False):
        
        # Model configuration
        self.model_name = "v_48_020"  # Which pretrained model to use
        self.backbone_noise = 0.000  # Std dev of Gaussian noise to add to backbone atoms
        self.max_length = 200000  # Maximum sequence length to process
        self.ca_only = ca_only  # Whether to use only CA atoms for scoring

        # device
        self.device = torch.device("cuda:0" if (torch.cuda.is_available()) else "cpu")

        # Set seed
        torch.manual_seed(seed)
        random.seed(seed)
        np.random.seed(seed)   
        
        # Model configuration
        hidden_dim = 128
        num_layers = 3 
    
        # Model weights path
        pwd = Path(__file__).resolve().parent
        if self.ca_only:
            checkpoint_path = str(pwd / f'ca_model_weights/{self.model_name}.pt')
        else:
            checkpoint_path = str(pwd / f'vanilla_model_weights/{self.model_name}.pt')

        # Load model
        checkpoint = torch.load(checkpoint_path, map_location=self.device) 
        self.model = ProteinMPNN(
            num_letters=21, 
            node_features=hidden_dim, 
            edge_features=hidden_dim, 
            hidden_dim=hidden_dim, 
            num_encoder_layers=num_layers, 
            num_decoder_layers=num_layers, 
            augment_eps=self.backbone_noise, 
            k_neighbors=checkpoint['num_edges'],
            ca_only=self.ca_only
        )
        self.model.to(self.device)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.model.eval()
        
    
    def _score_sequence(self, X, S, mask, chain_M, chain_M_pos, residue_idx, 
                    chain_encoding_all, sequence=None, num_batches=1,
                    decoding_order=None):

        global_scores = []
        alphabet_dict = dict(zip('ACDEFGHIKLMNPQRSTVWYX', range(21))) 
        
        # If FASTA sequence provided, update the input sequence tensor
        if sequence is not None:
            seq = sequence.replace('/', '')
            S[:,:len(seq)] = torch.tensor(
                [alphabet_dict[AA] for AA in seq], 
                device=self.device
            )[None,:].repeat(X.shape[0], 1)
        
        # Score sequences in batches
        for _ in range(num_batches):

            randn = torch.randn(chain_M.shape, device=X.device)
            log_probs = self.model(X, S, mask, chain_M*chain_M_pos, residue_idx, 
                                   chain_encoding_all, randn,
                                   decoding_order=decoding_order)

            # Calculate scores
            criterion = torch.nn.NLLLoss(reduction='none')
            loss = criterion(
                log_probs.contiguous().view(-1, log_probs.size(-1)),
                S.contiguous().view(-1)
            ).view(S.size())
            global_score = torch.sum(loss * mask, dim=-1) / torch.sum(mask, dim=-1)
            global_scores.append(global_score.cpu().numpy())

        score = np.array(global_scores)
        return score.mean(axis=0)

    
    def score(self,
            sequence: str | list[str] | None,
            structure: str | md.Trajectory,
            chains: str,
            num_seq_per_target: int = 1,
            relative: bool = True,
            decoding_order = None
        ) -> np.ndarray:

        all_chain_list = chains.split()

        if isinstance(structure, str):
            pdb_dict_list = parse_file(structure, all_chain_list, ca_only=self.ca_only) #[:10]
        else:
            pdb_dict_list = parse_traj(structure, all_chain_list, ca_only=self.ca_only)
        protein = StructureDatasetPDB(pdb_dict_list)

        # Process the input
        with torch.no_grad():
            # Prepare batch
            batch_clones = [p.copy() for p in protein]
            featurize_output = tied_featurize(batch_clones, self.device, all_chain_list, ca_only=self.ca_only)
            X, S, mask, chain_M, chain_encoding_all, chain_M_pos, residue_idx = featurize_output

            # score native seq only
            if sequence is None:
                pdb_score = self._score_sequence(
                    X, S, mask, chain_M, chain_M_pos, residue_idx, 
                    chain_encoding_all, sequence=None, 
                    num_batches=num_seq_per_target,
                    decoding_order=decoding_order
                )
                pdb_score = pdb_score[:, None]
                return pdb_score.T

            if relative:
                # Score native PDB sequence
                pdb_score = self._score_sequence(
                    X, S, mask, chain_M, chain_M_pos, residue_idx, 
                    chain_encoding_all, sequence=None, 
                    num_batches=num_seq_per_target,
                    decoding_order=decoding_order
                )
                pdb_score = pdb_score[:, None]
            else:
                pdb_score = np.zeros((X.shape[0], 1))
            
            if not isinstance(sequence, list):
                sequence = [sequence]
            
            query_score = np.array([
                self._score_sequence(
                    X, S, mask, chain_M, chain_M_pos, residue_idx,
                    chain_encoding_all, sequence=seq,
                    num_batches=num_seq_per_target,
                    decoding_order=decoding_order
                ) for seq in sequence
            ]).T

            delta_mean = (query_score - pdb_score) * len(sequence[0])

            return delta_mean.T
