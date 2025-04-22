import argparser
import synchronize

def run_pair_mode(args, params):
    print('Running SynchDP in pair mode...')
    from synchdpmod import synchDP
    synch_sev = synchronize.Pair(item, self.Rs, self.params.penalty, self.params.pcut, window, self.params.ecut)
    
    
def run_reference_mode(args, params):
    print('Running SynchDP in reference guided mode...')
    ref = synchronize.Reference(target = args.target,
                    query = args.query,
                    date_info = args.date_info,
                    omics_info = args.omics_info,
                    outdir = args.outdir,
                    params = params)


def run_denovo_mode(args, params):
    print('Running SynchDP in de novo mode...')
    net = synchronize.Denovo(query = args.query,
                  cover_ratio = args.cover_ratio,
                  date_info = args.date_info,
                  omics_info = args.omics_info,
                  outdir = args.outdir,
                  params = params)


if __name__ == '__main__':
    args, params=argparser.get()

    if args.mode == 'pair':
        run_pair_mode(args, params)

    if args.mode == 'reference':
        from synchronize import Reference
        if not args.target:
            raise ValueError('For reference guided mode, target sequence with \'-t\' option must be provided.')
        run_reference_mode(args, params)
        
    elif args.mode == 'denovo':
        from synchronize import Denovo
        run_denovo_mode(args, params)
        
    print('SynchDP complete.')
    