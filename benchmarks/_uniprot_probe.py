import urllib.request

JOB = "ncbiblast-R20261022-115340-0199-60127337-p1m"
url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/%s/json" % JOB
print(urllib.request.urlopen(url, timeout=60).read().decode()[:2500])
