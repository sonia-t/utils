# adapted from stinkpot
# does not work with some files with special chars
# should add check for whether ~/.Trash dir exists

for file in $*; do
    if [ -e $file ]; then

        # Target exists and can be moved to Trash safely
        if [ ! -e ~/.Trash/$file ]; then
            mv $file ~/.Trash

        # Target exists and conflicts with target in Trash
        elif [ -e ~/.Trash/$file ]; then

            # Increment target name until 
            # there is no longer a conflict
            i=1
            while [ -e ~/.Trash/$file.$i ];
            do
                i=$(($i + 1))
            done

            # Move to the Trash with non-conflicting name
            mv $file ~/.Trash/$file.$i
        fi

    # Target doesn't exist, return error
    else
        echo "rm: $file: No such file or directory";
    fi
done
