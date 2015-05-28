try:
	# Used BeautifulSoup 4.3.2 for HTML parsing and mechanize/cookielib for handling
	# login to GEDmatch
	from bs4 import BeautifulSoup
	import mechanize, cookielib, html2text

	# For building covariance matrix
	import numpy as np

	# For accepting command line arguments
	import argparse

	# For storing the covariance matrix
	import cPickle as pickle

	# For waiting when the website seems to be down temporarily
	import time

	# For exiting with error codes
	import sys
except Exception as e:
	print 'Error loading modules:\n%s' % e
	sys.exit(1)

def setup_browser():

	br = mechanize.Browser()
	cj = cookielib.LWPCookieJar()
	br.set_cookiejar(cj)
	br.set_handle_equiv(True)
	br.set_handle_redirect(True)
	br.set_handle_referer(True)
	br.set_handle_robots(False)
	br.set_handle_refresh(mechanize._http.HTTPRefreshProcessor(), max_time=1)
	br.add_headers = [('User-agent', 'Chrome')]
	return br

def login_to_gedmatch(br, email, password):

	try:
		br.open('http://v2.gedmatch.com/login1.php')
		br.select_form(nr=0)
		br.form['email'] = email
		br.form['password'] = password
		br.submit()
		return br
	except Exception as e:
		print 'Error while logging into GEDmatch:\n%s' % e
		print 'This is likely due to a temporary outage: will wait 10 s and reattempt.'
		return None

def login_to_gedmatch_wrapper(br, email, password):

	while True:
		br = login_to_gedmatch(br, email, password)
		if br != None:
			break
		time.sleep(10)
	return br

def get_first_n_hits(br, gedmatch_kit,  n):

	try:
		gedmatch_kits = []

		one_to_many_link = br.find_link(url='./r-list1z.php')
		br.follow_link(one_to_many_link)
		br.select_form(nr=0)
		br.form['kit_num'] = gedmatch_kit
		br.submit()

		soup = BeautifulSoup(br.response().read())
		matches = soup.find_all('tr')
		del matches[0:3] # These were headers

		for match in matches[0:n]:
			cells = match.find_all('td')
			gedmatch_kits.append(cells[0].contents[0].encode('ascii', 'ignore'))
		return gedmatch_kits
	except Exception as e:
		print 'Error while performing the one-to-many comparison:\n%s' % e
		print 'This is likely due to a temporary outage: will wait 10s and reattempt.'
		return None

def get_first_n_hits_wrapper(gedmatch_kit, n, email, password):

	br = setup_browser()

	while True:

		br = login_to_gedmatch_wrapper(br, email, password)
		gedmatch_kits = get_first_n_hits(br, gedmatch_kit,  n)
		if gedmatch_kits is not None:
			break
		time.sleep(10)

	return gedmatch_kits


def get_name_and_total_cMs(br, gedmatch_kit1, gedmatch_kit2):
	try:
		one_to_one_link = br.find_link(url='./u_compare1.php')
		br.follow_link(one_to_one_link)
		br.select_form(nr=0)
		br.form['kit1'] = gedmatch_kit1
		br.form['kit2'] = gedmatch_kit2
		br.submit()

		html = br.response().read()
		name = html.split(') and ', 1)[1].split(')<br>', 1)[0].split(' (')[1].encode('ascii', 'ignore')
		cMs = float(html.split('Total of segments > 7 cM = ')[1].split(' cM<br>')[0])	
		return cMs, name
	except Exception as e:
		print 'Error during one-to-one comparison of %s and %s:\n%s' % (gedmatch_kit1,
			gedmatch_kit2, e)
		print 'This is likely due to a temporary website outage: will wait 10s and reattempt.'
		return None, None

def get_name_and_total_cMs_wrapper(br, gedmatch_kit1, gedmatch_kit2, email, password):
	while True:
		br = login_to_gedmatch_wrapper(br, email, password)
		cMs, name = get_name_and_total_cMs(br, gedmatch_kit1, gedmatch_kit2)
		if cMs is not None:
			break
		time.sleep(10)
	return cMs, name

def build_cov_matrix(gedmatch_kits, email, password):
	n = len(gedmatch_kits)
	cov_matrix = np.zeros((n, n))
	names = ['']*n
	br = setup_browser()

	for i in xrange(0, len(gedmatch_kits)):
		print 'Currently running one-to-one comparisons for kit ' + \
			'%d of %d (%s)' % (i+1, n, gedmatch_kits[i])
		cov_matrix[i,i] = 1 # By definition

		for j in xrange(i+1, n):
			cMs, name = get_name_and_total_cMs_wrapper(br, gedmatch_kits[i],
				gedmatch_kits[j], email, password)
			names[j] = name # Will assign this many times; easier than the alternative
			fraction_autosomal_cMs_shared = cMs / (2 * 3578.1)
			cov_matrix[i,j] = fraction_autosomal_cMs_shared
			cov_matrix[j,i] = fraction_autosomal_cMs_shared

	names.pop(0)

	return cov_matrix, names

def get_segments(br, poi_kit, match_kit):

	try:
		segments = []

		# Browser is assumed to be on the select page (which comes right after
		# logging in)
		one_to_one_link = br.find_link(url='./u_compare1.php')
		br.follow_link(one_to_one_link)
		br.select_form(nr=0)
		br.form['kit1'] = poi_kit
		br.form['kit2'] = match_kit
		br.submit()

		# Now we parse the results. First we get the name corresponding to the
		# match kit.
		html = br.response().read()
		name = html.split(') and ', 1)[1].split(')<br>', 1)[0].split(' (')[1].encode('ascii', 'ignore')

		# Now we get the info on each segment
		soup = BeautifulSoup(html)
		rows = soup.find_all('tr')
		rows.pop(0) # This was just a header
		for row in rows:
			cells = row.find_all('td')
			chromosome = cells[0].contents[0].encode('ascii', 'ignore')
			start_pos = int(cells[1].contents[0])
			end_pos = int(cells[2].contents[0])
			cMs = float(cells[3].contents[0])
			SNPs = int(cells[4].contents[0])
			segments.append([chromosome, start_pos, end_pos, cMs, SNPs, match_kit, name])

		return segments

	except Exception as e:
		print 'Error while comparing %s to %s:\n%s' % (poi_kit, match_kit, e)
		print 'This is likely due to a temporary outage: will wait 10s and reattempt.'
		return None

def get_segments_wrapper(br, poi_kit, match_kit, email, password):
	while True:
		br = login_to_gedmatch_wrapper(br, email, password)
		segments = get_segments(br, poi_kit, match_kit)
		if segments is not None:
			break
		time.sleep(10)
	return segments

def get_segments_multiple_ids(poi_kit, gedmatch_kits, email, password):

	# Start by finding all segments shared between the person of interest and
	# all of the other gedmatch kits.
	segments = []
	br = setup_browser()
	for match_kit in gedmatch_kits:
		results = get_segments_wrapper(br, poi_kit, match_kit, email, password)
		if len(results) > 0:
			segments.extend(results)

	segments.sort(key=lambda x: (x[0], x[1], x[2]))

	return segments

def main(email, password, gedmatch_kit, n, output_filename):

	print 'Beginning with a one-to-many comparison of GEDmatch kit ' + \
		'%s against all others to find the top %d best hits. This is slow.' % (gedmatch_kit,
		n)
	gedmatch_kits = get_first_n_hits_wrapper(gedmatch_kit, n, email, password)
	all_kits = [gedmatch_kit]
	all_kits.extend(gedmatch_kits)

	print 'Having identified the top hits, proceeding to covariance matrix construction.'
	cov_matrix, names = build_cov_matrix(all_kits, email, password)

	print '\n*** Printing the covariance matrix ***'
	print cov_matrix

	print '\n*** Printing the GEDmatch kit ids and names in the same order as ' + \
		'the rows and columns of the covariance matrix ***'
	print 'GEDmatch kit %s: proband' % gedmatch_kit
	for i in xrange(0,len(names)):
		print 'GEDmatch kit %s: %s' %(gedmatch_kits[i], names[i])

	print '\nCompiling a list of segments shared with the top hits'
	segments = get_segments_multiple_ids(gedmatch_kit, gedmatch_kits, email, password)
	print '\n*** Printing the segment data. Note that the total number of ' + \
		'cMs for these segments may be less than the number of cMs reported ' + \
		'as shared in the information used to construct the covariance matrix ***'
	for segment in segments:
		print 'Chr %s:%d-%d (%0.1f cM/%d SNPs) matches %s (%s)' % tuple(segment)

	try:
		output_file = open(output_filename, 'wb')
		pickle.dump(cov_matrix, output_file)
		pickle.dump(names, output_file)
		pickle.dump(segments, output_file)
	finally:
		output_file.close()

	return

if __name__ == "__main__":

	'''
	Example usage:

	python build_cov_matrix_from_gedmatch.py -e <your email> -p <your password>
		-n 3 -g M680541 -o test.pickle

	Note that you can use any account to log in for this purpose; you do not
	need to use the one that uploaded the GEDmatch kit.
	'''

	parser = argparse.ArgumentParser(description='Finds the top n best matches ' + \
		'for an input GEDmatch kit number. Performs all pairwise comparisons to ' + \
		'estimate the similarity of each match to the proband and each other. ' + \
		'Outputs a similarity matrix (autosomal cMs shared/total autosomal cMs) ' + \
		'as well as the names and GEDmatch ids of the rows and columns.')
	parser.add_argument('--email', '-e', type=str, required=True,
		help='email address of GEDmatch account', dest='email', action='store')
	parser.add_argument('--password', '-p', type=str, required=True,
		help='password of GEDmatch account', dest='password', action='store')
	parser.add_argument('--number_of_matches', '-n', type=int, required=True,
		help='number of top matches to consider', dest='n', action='store')
	parser.add_argument('--gedmatch_id', '-g', type=str, required=True,
		help='GEDmatch kit number', dest='gedmatch_kit')
	parser.add_argument('--output_filename', '-o', type=str, required=True,
		help='Filename for storing covariance matrix (in pickled form)',
		dest='output_filename', action='store')

	try:
		args = parser.parse_args()
		main(args.email, args.password, args.gedmatch_kit, args.n, args.output_filename)
	except Exception as e:
		print 'GEDmatch parsing failed:\n%s' % e
		sys.exit(2)
