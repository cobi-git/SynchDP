from dataclasses import dataclass
import argparse

@dataclass
class Params:
    penalty: float
    pcut: float
    window: int
    ecut: float

def get():
	# mandatory parameters
	parser = argparse.ArgumentParser(description='Choose between pair, de novo and reference guided mode.')
	parser.add_argument('-m', '--mode', choices=['pair', 'denovo', 'reference'], required=True)
	parser.add_argument('-q', '--query', type=str, required=True, help='The input sequence(s) to align.')

	# optional parameters
	parser.add_argument('-t', '--target', type=str, help='The target sequence (or reference sequence).')
	parser.add_argument('-r','--cover_ratio', type=int, default=95, help='For denovo mode only.')
	parser.add_argument('--date_info', type=str, default=None, help='Optional: clinical date information.')
	parser.add_argument('--omics_info', type=str, default=None, help='Optional: omics sampling information.')
	parser.add_argument('-o', '--outdir', type=str, default='./', help='The output path of the results')

	# common parameters
	parser.add_argument('-p', '--penalty', type=float, default=0.01, help='The penalty weight for time gap [float: 0~1]')
	parser.add_argument('-c', '--pcut', type=float, default=0.7, help='The minimum required correlation of the alignment [float: 0~1]')
	parser.add_argument('-w', '--window', type=int, default=6, help='The length of window used for alignment [int]')
	parser.add_argument('-e', '--ecut', type=float, default=3.5, help='The minimum required Euclidean distance of the alignment. May need to be decided based on the input data. [float]')

	args = parser.parse_args()

	params = Params(
	    penalty=args.penalty,
	    pcut=args.pcut,
	    window=args.window,
	    ecut=args.ecut)

	return args, params