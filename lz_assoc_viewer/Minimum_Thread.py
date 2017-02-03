import threading
import minimum_solution
exitFlag = 0

class Minimum_Thread(threading.Thread):
	def __init__(self, filename, filetype, numMin, minDist):
		threading.Thread.__init__(self)
		self.filename = filename
		self.filetype = filetype
		self.numMin = numMin
		self.minDist = minDist
		self.minimums = None
	def run(self):
		print("starting minimum thread....")
		self.minimums = minimum_solution.find_min_pvals(self.filename, self.filetype, self.numMin, self.minDist)
		print("finishing minimum thread...")
	def join(self):
		threading.Thread.join(self)
		return self.minimums





