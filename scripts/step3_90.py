# STEP 3 STARTS HERE:
# make a non-redundant list from the sets
# generated in step 2

import json

final_list = []
with open('bins2/parallel_merged.out' ,'r') as in_file:
	master_list = json.load(in_file)

while True:
	test = master_list[0]
	working_list = test
	for test_list in master_list[1:]:
		if not (set(test_list).intersection(test)):
			a = False
		else:
			a = True
		if a == False:
			next
		else:
			working_list.extend(test_list)
			master_list.remove(test_list)
	master_list.remove(test)
	x = set(working_list)
	x = list(x)
	final_list.append(x)
	if master_list == []:
		break

with open('bins2/total_step3_list.out', 'w') as out_file:
	json.dump(final_list, out_file)

