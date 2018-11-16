#=============================
# read in command line args
# space separated
#=============================
while [[ $# > 1 ]]
do
    key="$1"
    shift
    case $key in
	--rnaseq_fq)
	    rnaseq_fq="$1"
	    shift
	    ;;
	--riboseq_fq)
	    riboseq_fq="$1"
	    shift
	    ;;
	--transcript_fa)
	    transcript_fa="$1"
	    shift
	    ;;
	--contaminant_fa)
	    contaminant_fa="$1"
	    shift
	    ;;
	--cds_range)
	    cds_range="$1"
	    shift
	    ;;
	--work_dir)
	    work_dir="$1"
	    # star index
	    star_idx_dir=${work_dir}StarIndex/
	    # star outputs
	    tmp_dir=${work_dir}alignment/
	    # salmon
	    sm_odir=${work_dir}sm_quant
	    # ribomap
	    output_dir=${work_dir}outputs
	    shift
	    ;;
	--nproc)
	    nproc="$1"
	    shift
	    ;;
	--adapter)
	    adapter="$1"
	    shift
	    ;;
	--min_fplen)
	    min_fplen="$1"
	    shift
	    ;;
	--max_fplen)
	    max_fplen="$1"
	    shift
	    ;;
	--offset)
	    offset="$1"
	    shift
	    ;;
	--nmismatch)
	    nmismatch="$1"
	    shift
	    ;;
	--tabd_cutoff)
	    tabd_cutoff="$1"
	    shift
	    ;;
	--star_idx_dir)
	    star_idx_dir="$1"
	    shift
	    ;;
	--alignment_dir)
	    tmp_dir="$1"
	    shift
	    ;;
	--sailfish_dir)
	    sm_odir="$1"
	    shift
	    ;;
	--output_dir)
	    output_dir="$1"
	    shift
	    ;;
	--force)
	    force="$1"
	    shift
	    ;;
	--useSecondary)
	    useSecondary="$1"
	    shift
	    ;;
	--rnaUnstranded)
	    useRC=true
	    shift
	    ;;
	--softClipping)
	    softclip="$1"
	    shift
	    ;;
	*)
            # unknown option
	    echo "unrecognized option $key!"
	    ;;
    esac
done