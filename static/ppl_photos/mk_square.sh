
# Get width and height
read w h < <(identify -format "%w %h" adrian.png)

# Check them
echo $w,$h

# Set `n` to lesser of width and height
n=$w
[ $h -lt $n ] && n=$h

# Now do actual crop
convert adrian.png -gravity center -extent "${n}x${n}" result.jpg
