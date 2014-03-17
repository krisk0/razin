#!/bin/sh

# Vile magic to strip nmod_mat_HNF-debug.c of debug code and repeated empty 
#  lines

m='#define NDEBUG 1'
m="This C code is automagically generated out of nmod_mat_HNF-debug.c\n\n$m"

perl -pe 's/ *#/#/g' nmod_mat_HNF-debug.c |
perl -pe s/\\x0a/\\x1a/g |
perl -pe 's/#if .*?#endif\x1a//g' |
perl -pe s/\\x1a\\x1a\\x1astatic/\\x1a\\x1astatic/g |
perl -pe s/\\x1a/\\x0a/g | 
sed -e "s=#define BUG.*=//$m=" | sed -e 's:#define BUG.*::' |
sed -e /MMUL/d | sed -e /MPLUS/d | sed -e '/delete me/d' \
> nmod_mat_HNF.c

#perl -pe s/\\x0a\\x0astatic/=static/g -i nmod_mat_HNF.c

cat << EOF >> nmod_mat_HNF.c
/*
I will be greatly irritated if you tell me that this algorithm is wrong
 WITHOUT GIVING EXAMPLE OF INPUT DATA that make it fail

Report bugs via Github mechanism or e-mail

My e-mail is in my blog, detailed information on how to get it is close to tail
 of setup.py
*/
EOF

# perl -pe s/\\x0a/\\x1a/g | 
# perl -p -e s:\\x1a\\x1a\\x1a\\x1a:\\x1a\\x1a\\x1a: | perl -p -e s:\\x1a\\x1a:\\x1a:g |
# perl -pe s/\\x1a/\\x0a/g \
