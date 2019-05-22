sed -n -r "s/.*integral = ([0-9.-]+).*/\1/p" output.txt

# if [ -f "data.txt" ]; then do
#   echo data.txt already exists
#   exit 
# done
#
# for x in {0..3}; do
#   file=x$x.txt
#   echo "reading from $file"
#   sed -n -r "s/.*integral = ([0-9.-]+).*/\1/p" $file >> data.txt
# done
#
# echo "number of lines in data.txt $(cat data.txt | wc -l)"

