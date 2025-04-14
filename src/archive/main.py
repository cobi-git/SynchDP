import argparse
import pickle
import pandas as pd
from dataclasses import dataclass
from denovo_module import DenovoAnalysis
from reference_module import ReferenceAnalysis


@dataclass
class Params:
    penalty: float
    pcut: float
    window: int
    ecut: float


def run_denovo_mode(querry_set, cover_ratio, date_info=None, omics_info=None, outdir=None, params = None):
    print("Running SynchDP in de novo mode...")
    net = DenovoAnalysis(querry_set = querry_set,
                          cover_ratio = cover_ratio,
                          date_info = date_info,
                          omics_info = omics_info,
                          outdir = outdir,
                          params = params)
    print("SynchDP complete.")

    
def run_reference_mode(reference_seq, querry_set, date_info=None, omics_info=None, outdir=None, params = None):
    print("Running SynchDP in reference guided mode...")
    ref = ReferenceAnalysis(reference_seq = reference_seq,
                            querry_set = querry_set,
                            date_info = date_info,
                            omics_info = omics_info,
                            outdir = outdir,
                            params = params)
    print("SynchDP complete.")

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Choose between de novo and reference guided mode.")
    parser.add_argument("--mode", choices=["denovo", "reference"], required=True)
    parser.add_argument("--querry_set", type=str, required=True)
    parser.add_argument("--outdir", type=str, default='./', help='The output path of the results')

    # 선택적 인자들
    parser.add_argument("--cover_ratio", type=int, default=95, help="For denovo mode only.")
    parser.add_argument("--reference_seq", type=str, help="For reference mode only.")
    parser.add_argument("--date_info", type=str, default=None, help="Optional: clinical date information.")
    parser.add_argument("--omics_info", type=str, default=None, help="Optional: omics sampling information.")

    # Shared SynchDP hyperparameters
    parser.add_argument("--penalty", type=float, default=0.01)
    parser.add_argument("--pcut", type=float, default=0.7)
    parser.add_argument("--window", type=int, default=6)
    parser.add_argument("--ecut", type=float, default=3.5)
    
    args = parser.parse_args()
    
    params = Params(
        penalty=args.penalty,
        pcut=args.pcut,
        window=args.window,
        ecut=args.ecut)
    
    if args.mode == "reference":
        if not args.reference_seq:
            raise ValueError("For reference mode, --reference_seq must be provided.")
        run_reference_mode(args.reference_seq, args.querry_set, args.date_info, args.omics_info, args.outdir, params)
        
    elif args.mode == "denovo":
        run_denovo_mode(args.querry_set, args.cover_ratio, args.date_info, args.omics_info, args.outdir, params)
        
        
    