#!/bin/bash
set -e -o pipefail
# Fetch the shared tools into `deps/` and record them in `deps/tools_version`. They are fetched from the
# JuliaPackageTools repository on GitHub, or from another clone of it given as the first argument.
#
# This script is one of the tools, so running it from `deps/` replaces it while it runs. Bash reads a script as it
# goes, so the whole of it is a single block, which bash reads before running any of it, and which exits at its end.
{
    SOURCE=${1:-https://github.com/tanaylab/JuliaPackageTools.git}
    TEMPORARY=`mktemp -d`
    trap "rm -rf $TEMPORARY" EXIT
    git clone --quiet --depth 1 $SOURCE $TEMPORARY 2> >(grep -v 'depth is ignored in local clones' >&2)
    FILES=`git -C $TEMPORARY ls-files | grep -v '^README.md$'`
    if [ -f deps/tools_version ]
    then
        for FILE in `tail -n +2 deps/tools_version | sed 's/^[^ ]*  //'`
        do
            if ! echo "$FILES" | grep -q -x "$FILE"
            then
                # The deletion is staged, since the Makefile lists the files by `git ls-files`, which would still
                # list a file only deleted from the disk.
                git rm -q -f --ignore-unmatch deps/$FILE
                rm -f deps/$FILE
            fi
        done
    fi
    for FILE in $FILES
    do
        mkdir -p deps/`dirname $FILE`
        cp -p $TEMPORARY/$FILE deps/$FILE
    done
    (echo "commit `git -C $TEMPORARY rev-parse HEAD`"; cd deps && sha256sum $FILES) > deps/tools_version
    exit
}
