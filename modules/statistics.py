#/usr/bin/python

import sqlite3
import time
import subprocess


#db = sys.argv[1]
#ch = "1"
#pos = 17900000

def get_stats(db, ch, pos1, pos2, typ):

	start = time.time()

	ch = ch.upper()
	ch = '{}{}'.format('chr', ch)
	print 'Normalizing the read coverage and generating statistics for this breakpoint'

	# create/define dictionary to save all stats in 
	statistics = {}

	# Connect the two databases
	db1 = sqlite3.connect(db)
	cursor = db1.cursor()


	# Classify SV using statistics for breakpoints:
	pos1 = int(pos1)/100*100
	pos2 = int(pos2)/100*100 

	#========================================================
	#  Get region specific information from SQLite
	#========================================================

	def region_info(ch, pos1, pos2):
			# If positions in reversed order, reverse it back
		if pos1 > pos2:
			a = pos1
			b = pos2
			pos1 = b
			pos2 = a

		print ch, pos1, pos2 	
		# get out the specific statistic.  
		cursor.execute('''SELECT avg(mappability_score) FROM table1 where chrom=? and start_pos between ? and ?''', (ch, pos1, pos2))
		map_i = cursor.fetchone()
		if map_i[0] is None:
			print 'map_i is None', map_i
			return 'NA', 'NA', 'NA'
		else:
			map_i = map_i[0]	

		cursor.execute('''SELECT GC_content FROM table1 where chrom=? and start_pos = ?''', (ch, pos1))
		gc_content = cursor.fetchone()
		if gc_content is None:
			print 'gc_content is None', gc_content
			return 'NA', 'NA', 'NA'
		else:
			gc_content = gc_content[0]	

		cursor.execute('''SELECT avg(read_coverage) FROM table1 where chrom=? and start_pos between ? and ?''', (ch, pos1, pos2))
		r_i = cursor.fetchone()
		if r_i[0] is None:
			print 'RD_i is None', r_i
			return 'NA', 'NA', 'NA'
		else:
			r_i = r_i[0]	

		return map_i, gc_content, r_i	


	RD_all = cursor.execute('''SELECT avg(read_coverage) FROM table1 where chrom = ? and read_coverage > ? and mappability_score > ?''', (ch, 0, 0.25)).fetchone()[0]

	
	t1 = time.time() - start


	#===================================================================
	#  Normalize read depth for all three regions by getting region info and calculating s
	#===================================================================
	start2 = time.time()
	if typ == 'same':
		print ch, pos1, pos2
		# Region 1
		map_1, gc_1, RD_1 = region_info(ch, pos1, pos2)
		print map_1, gc_1, RD_1
		if map_1 == 'NA':
			return statistics
		RD_gc_1 = cursor.execute('''SELECT avg(read_coverage) FROM table1 where chrom = ? and GC_content=? and read_coverage > ?''', (ch, gc_1, 0)).fetchone()[0]
		RD_norm_1 = RD_1 * RD_all / RD_gc_1
		t2 = time.time() - start2

		start3 = time.time()
		# Region 2 (check the region around breakpoint 1)
		bp1_1 = pos1 - 500
		bp1_2 = pos1 + 500
		map_2, gc_2, RD_2 = region_info(ch, bp1_1, bp1_2)
		if map_2 == 'NA':
			return statistics	
		RD_gc_2 = cursor.execute('''SELECT avg(read_coverage) FROM table1 where chrom = ? and GC_content=? and read_coverage > ?''', (ch, gc_2, 0)).fetchone()[0]
		RD_norm_2 = RD_2 * RD_all / RD_gc_2
		t3 = time.time() - start3

		start4 = time.time()
		# Region 3 (check the region around breakpoint 2)
		bp2_1 = pos2 - 500
		bp2_2 = pos2 + 500
		map_3, gc_3, RD_3 = region_info(ch, bp2_1, bp2_2)
		if map_3 == 'NA':
			return statistics	
		RD_gc_3 = cursor.execute('''SELECT avg(read_coverage) FROM table1 where chrom = ? and GC_content=? and read_coverage > ?''', (ch, gc_3, 0)).fetchone()[0]

		RD_norm_3 = RD_3 * RD_all / RD_gc_3
		t4 = time.time() - start4

	elif typ != 'same':
		t2 = 0
		start3 = time.time()
		# Region 2 (check the region around breakpoint 1)
		bp1_1 = pos1 - 500
		bp1_2 = pos1 + 500
		map_2, gc_2, RD_2 = region_info(ch, bp1_1, bp1_2)
		if map_2 == 'NA':
			return statistics	
		RD_gc_2 = cursor.execute('''SELECT avg(read_coverage) FROM table1 where chrom = ? and GC_content=? and read_coverage > ?''', (ch, gc_2, 0)).fetchone()[0]
		RD_norm_2 = RD_2 * RD_all / RD_gc_2
		t3 = time.time() - start3

		start4 = time.time()
		# Region 3 (check the region around breakpoint 2)
		bp2_1 = pos2 - 500
		bp2_2 = pos2 + 500
		chrom2 = typ
		map_3, gc_3, RD_3 = region_info(chrom2, bp2_1, bp2_2)
		if map_3 == 'NA':
			return statistics	
		RD_gc_3 = cursor.execute('''SELECT avg(read_coverage) FROM table1 where chrom = ? and GC_content=? and read_coverage > ?''', (ch, gc_3, 0)).fetchone()[0]

		RD_norm_3 = RD_3 * RD_all / RD_gc_3
		t4 = time.time() - start4		
			

	subprocess.call('echo '+ str(t1) + ' ' + str(t2) + ' ' + str(t3) + ' ' + str(t4) + ' >> time_1.txt', shell = True)

	# Save all info in the dictionary and return for SV classificaton and genotyping
	if typ == 'same':
		statistics['RD_1'] = RD_1 # read count for the region
		statistics['RD_gc_1'] = RD_gc_1 # average read count of all windows having the same gc percentage as i:th window
		statistics['gc_1'] = gc_1 # gc content (refernece genome) for the i:th window
		statistics['RD_norm_1'] = RD_norm_1 # normalized read counts for the i:th window
		statistics['map_1'] = map_1 # mappability score for the ith window


	statistics['RD_2'] = RD_2 # read count for the region
	statistics['RD_3'] = RD_3 # read count for the region
	statistics['RD_all'] = RD_all # average read count of all windows with a value higher then 0 
	statistics['RD_gc_2'] = RD_gc_2 # average read count of all windows having the same gc percentage as i:th window
	statistics['RD_gc_3'] = RD_gc_3 # average read count of all windows having the same gc percentage as i:th window
	statistics['gc_2'] = gc_2 # gc content (refernece genome) for the i:th window
	statistics['gc_3'] = gc_3 # gc content (refernece genome) for the i:th window
	statistics['RD_norm_2'] = RD_norm_2 # normalized read counts for the i:th window
	statistics['RD_norm_3'] = RD_norm_3 # normalized read counts for the i:th window
	statistics['map_2'] = map_2 # mappability score for the ith window
	statistics['map_3'] = map_3 # mappability score for the ith window

	db1.commit()

	print 'found statistics'
	return statistics





