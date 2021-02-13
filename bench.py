import subprocess

binary = './findprimes'
test_count = 4
files = [ '8_test_half.txt' ]
thread_counts = [1, 2, 3, 4, 6, 8, 12]
block_size_mins = [1, 2, 4, 8]

tests = [ {'file': f, 'thread_count': c, 'block_size_min': b, 'durations': []}
	for f in files 
	for c in thread_counts
	for b in block_size_mins
]
for i in range(test_count):
	for test in tests:
		test['durations'] += [float(subprocess.run(
				[binary, str(test['thread_count']), test['file']],
				stderr=subprocess.PIPE
			).stderr)
		]

print('\n'.join([str(test) for test in tests]))
