#=============================
# read in command line args
# space separated
#=============================

while [[ $# > 1 ]]
do
    echo "number of args: $#"
    key="$1"
    shift
    case $key in
    --test)
        test="$1"
        shift
        ;;
    --test2)
        test2="$1"
        shift
        ;;
    *)
            # unknown option
        echo "unrecognized option $key!"
        ;;
    esac
done

echo "test: $test"
echo "test2: $test2"
