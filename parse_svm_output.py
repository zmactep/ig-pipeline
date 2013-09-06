__author__ = 'Kos'

with open('data.txt') as f:
    current_name = prediction = ''
    for line in f:
        name, val = line.rstrip().split('\t')[:2]
        if current_name == name:
            prediction += val
        else:
            print(prediction)
            current_name = name
            prediction = current_name + '\t' + val
