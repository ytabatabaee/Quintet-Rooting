with open('test_indices_with_lengths.txt', 'r') as fp:
    indices_list_str = fp.read()

indices_list_str = indices_list_str.split(')')
indices_list = []

for item in indices_list_str:
    indices_list.append([int(s) for s in item.replace('(', ' ').replace(',', ' ').split() if s.isdigit()])
indices_list = indices_list[:-1]
print(len(indices_list))


#print(indices_list_str)
