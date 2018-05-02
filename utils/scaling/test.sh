
function change_cfd() {
	sed -i 4s/.*/$1/ input
	sed -i 5s/.*/$2/ input
	sed -i 6s/.*/$3/ input
}

function change_md() {
	sed -i 7s/.*/$1/ input
	sed -i 8s/.*/$2/ input
	sed -i 9s/.*/$3/ input
}

function factors() {
	for i in $(seq 1 $1)
	do
	 [ $(expr $1 / $i \* $i) == $1 ] && echo $i
	done
}

MAX=96
for ((i=12;i<=$MAX;i=i*2)); do 
	for ((j=$i;j<=$MAX;j=j*2)); do
	    factors $i 
	    factors $j
	    echo aprun -n $i ./cfd : -n $j ./md
	done
done
#change_cfd 8 3 2
#change_md 8 3 2
