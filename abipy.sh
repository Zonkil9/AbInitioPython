echo "Do you want to save results into a file?"
read -p "Type \"yes\" or \"no\": " word
if [ "$word" = "yes" ]; then
    read -p "Provide filename of an output: " output
    script -c 'python main.py' $output
    sed -i '1d' $output
    sed -i '$d' $output
elif [ "$word" = "no" ]; then
    python main.py
else
    echo "You need to type \"yes\" or \"no\"!"
fi
