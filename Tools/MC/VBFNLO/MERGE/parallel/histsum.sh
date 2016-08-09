#! /bin/sh

# arguments: <output-dir> <input-dirs> 

if [ $# -le 1 ]; then
  echo "Usage: $0 <output-dir> <hist-dir> <log-file> <input-dirs> <...>" >&2
  exit 1
fi

outdir=$1
histdir=$2
logfile=$3

echo $outdir
echo $histdir
echo $logfile

shift 3

mkdir -p "$outdir/$histdir"
mkdir -p "$outdir/$histdir/LO"
mkdir -p "$outdir/$histdir/NLO"
mkdir -p "$outdir/$histdir/Kfac"

for i in $1/$histdir/LO/hist.*; do
histname=`basename $i`
histexe=`dirname $0`/histsum 
histout=`$histexe $histname $outdir $histdir "$@" 2>&1 `
echo "$histout"
histout=`echo "$histout" | grep "Cannot open file" | awk '{print $4}' | awk -F'/' '{print $1}' | sort | uniq`
done

for i in $1/$histdir/LO/hist2.*; do
histname=`basename $i`
histexe2=`dirname $0`/histsum2d 
histout=`$histexe2 $histname $outdir $histdir "$@" 2>&1 `
echo "$histout"
histout=`echo "$histout" | grep "Cannot open file" | awk '{print $4}' | awk -F'/' '{print $1}' | sort | uniq`
done

k=0
for i in $@; do
  echo $histout | grep $i >/dev/null 2>&1 && continue
  k=$((k+1))
  lineborn=`grep 'TOTAL result (LO):' $i/$logfile`
  linevirt1=`grep 'NLO virtual result:' $i/$logfile`
  linereal=`grep 'QCD real emission:' $i/$logfile`
  linetotl=`grep 'TOTAL result (NLO):' $i/$logfile`
  born=`echo "$lineborn" | awk '{print $4}'`
  real=`echo "$linereal" | awk '{print $4}'`
  totl=`echo "$linetotl" | awk '{print $4}'`
  virt1=`echo "$linevirt1" | awk '{print $4 "+"}' | tr -d '\n' `
  errborn=`echo "$lineborn" | awk '{print $6}'`
  errreal=`echo "$linereal" | awk '{print $6}'`
  errtotl=`echo "$linetotl" | awk '{print $6}'`
  errvirt1=`echo "$linevirt1" | awk '{print "(" $6 ")^2+"}' | tr -d '\n' `
  totborn="$totborn$born+"
  totreal="$totreal$real+"
  tottotl="$tottotl$totl+"
  totvirt1="$totvirt1$virt1"
  toterrborn="$toterrborn($errborn)^2+"
  toterrreal="$toterrreal($errreal)^2+"
  toterrtotl="$toterrtotl($errtotl)^2+"
  toterrvirt1="$toterrvirt1$errvirt1"
done
echo "Got $k entries from log files."

virt1=`echo "scale=10;${virt1}0" | sed 's.E.\\*10\\^.g' | bc | tr -d '~' `
errvirt1=`echo "scale=10;sqrt(${errvirt1}0)" | sed 's.E.\\*10\\^.g' | bc | tr -d '~' `

totborn=`echo "scale=10;(${totborn}0)/$k" | sed 's.E.\\*10\\^.g' | bc | tr -d '~' `
totreal=`echo "scale=10;(${totreal}0)/$k" | sed 's.E.\\*10\\^.g' | bc | tr -d '~' `
tottotl=`echo "scale=10;(${tottotl}0)/$k" | sed 's.E.\\*10\\^.g' | bc | tr -d '~' `
totvirt1=`echo "scale=10;(${totvirt1}0)/$k" | sed 's.E.\\*10\\^.g' | bc | tr -d '~' `
toterrborn=`echo "scale=10;sqrt(${toterrborn}0)/$k" | sed 's.E.\\*10\\^.g' | bc | tr -d '~' `
toterrreal=`echo "scale=10;sqrt(${toterrreal}0)/$k" | sed 's.E.\\*10\\^.g' | bc | tr -d '~' `
toterrtotl=`echo "scale=10;sqrt(${toterrtotl}0)/$k" | sed 's.E.\\*10\\^.g' | bc | tr -d '~' `
toterrvirt1=`echo "scale=10;sqrt(${toterrvirt1}0)/$k" | sed 's.E.\\*10\\^.g' | bc | tr -d '~' `
totratioborn=`echo "scale=10;${toterrborn}/${totborn}*100" | sed 's.E.\\*10\\^.g' | bc | tr -d '~' `
totratioreal=`echo "scale=10;${toterrreal}/${totreal}*100" | sed 's.E.\\*10\\^.g' | bc | tr -d '~' `
totratiototl=`echo "scale=10;${toterrtotl}/${tottotl}*100" | sed 's.E.\\*10\\^.g' | bc | tr -d '~' `
totratiovirt1=`echo "scale=10;${toterrvirt1}/${totvirt1}*100" | sed 's.E.\\*10\\^.g' | bc | tr -d '~' `
totkfac=`echo "scale=10;${tottotl}/${totborn}" | sed 's.E.\\*10\\^.g' | bc | tr -d '~' `

echo "TOTAL result (LO):   $totborn  +-  $toterrborn fb  $totratioborn % " >$outdir/LOG.out
echo "NLO virtual result:  $totvirt1  +-  $toterrvirt1 fb  $totratiovirt1 % " >>$outdir/LOG.out
echo "QCD real emission:   $totreal  +-  $toterrreal fb  $totratioreal % " >>$outdir/LOG.out
echo "TOTAL result (NLO):  $tottotl  +-  $toterrtotl fb  $totratiototl % " >>$outdir/LOG.out
echo "K-factor:            $totkfac " >>$outdir/LOG.out

