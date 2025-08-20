export LD_LIBRARY_PATH=$PWD/aces2:$LD_LIBRARY_PATH

echo "Do you want to save results into a file?"
read -r -p "Type \"yes\" or \"no\": " word
if [ "$word" = "yes" ]; then
    read -r -p "Provide filename of an output: " output
    echo "Preparing ACES2 files"
    mkdir jobscr
    rm jobscr/*
    cd jobscr || exit
    cp ../ZMAT ZMAT
    ln -s ../aces2/GENBAS GENBAS
    {
      ../aces2/xjoda
      ../aces2/xvmolmod
      ../aces2/xvprops
    } >> scratch.out
    cp ../python/main.py ../python/quantum.py .
    script -c 'python main.py' "$output"
    sed -i '1d' "$output"
    sed -i '$d' "$output"
    mv "$output" ..
    cd ..
    rm -rf jobscr
elif [ "$word" = "no" ]; then
    echo "Preparing ACES2 files"
    mkdir jobscr
    rm jobscr/*
    cd jobscr || exit
    cp ../ZMAT ZMAT
    ln -s ../aces2/GENBAS GENBAS
    {
      ../aces2/xjoda
      ../aces2/xvmolmod
      ../aces2/xvprops
    } >> scratch.out
    cp ../python/main.py ../python/quantum.py .
    python main.py
    cd ..
    rm -rf jobscr
else
    echo "You need to type \"yes\" or \"no\"!"
fi
