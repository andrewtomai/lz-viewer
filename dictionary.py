##working with dictionaries
import json
variant = "variant"
position = "position"

data = {variant : [], position : []}
data[variant].append("dubstep")
data[variant].append("wubstep")
data[variant].append("trubstep")
data[position].append(1)
data[position].append(2)
data[position].append(3)

json_data = json.dumps(data)
print json_data
