#!/bin/bash
IN=$(basename $1)
NAME=${IN%.cl}

echo "const char *"${NAME}" =" > $1.h
/usr/bin/sed -e 's/\\/\\\\/g;s/"/\\"/g;s/^/"/;s/$/\\n"/;s/real/double/g' \
    $1 >> $1.h
echo ";" >>$1.h
