mkdir -p bin
ls *.c | while read F; do gcc -o `echo $F | cut -d . -f 1 - ` $F -lcgns -lhdf5 -lm; mv `echo $F | cut -d . -f 1 - ` bin/; done

