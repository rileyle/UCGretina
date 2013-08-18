rm -f awalls
rm -f aclust
rm -f asolid
rm -f aslice
rm -f aeuler


ln -s ${1}clust.list aclust
ln -s ${1}walls.list awalls
ln -s ${1}solid.list asolid
ln -s ${1}slice.list aslice
ln -s ${1}euler.list aeuler

echo "geometry reconfigured to ${1}"

ls -lht a*
