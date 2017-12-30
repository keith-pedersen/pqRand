# This script validates the distributions of pYqRand, outputting the results to validation_report.py

import pYqRand as pqr
import numpy
import scipy.integrate as integ
import math

def Banner(text, char = "="):
	banner = char*len(text) + "\n"
	banner += text + "\n"
	banner += char*len(text) + "\n\n"
	return banner
	
def BinaryAccumulate_Destructive(array):
	if(len(array) % 2):
		array = numpy.append(array, 0.)
	
	size = len(array)
	
	while(size > 1):
		if(size % 2):
			array[size] = 0.
			size += 1
		size //= 2
		for i in range(size):
			array[i] += array[i + size]
	
	return array[0]
	
def Squared(x):
	return x*x
	
def BinaryAccumulate(array):
	return BinaryAccumulate_Destructive(numpy.copy(array))

def TestDistribution(dist, sampleSize, binSize):
	report = ""
	
	report += Banner("Validation Report: " + str(dist), char = "#")
	
	gen = pqr.engine()
	indexer = pqr.uniform_integer(0, sampleSize)
		
	sample = dist.GetSample(sampleSize, gen)
	
	assert(len(sample) == int(sampleSize))
	mean = BinaryAccumulate(sample) / sampleSize
	
	squareDeviations = numpy.asarray(list(map(lambda x: Squared(x - mean), sample)))
	variance = BinaryAccumulate_Destructive(squareDeviations) / (sampleSize - 1.)
	
	#~ sample = numpy.insert(sample, 0, [dist.min(), dist.max()])
	sample = numpy.sort(sample)
	
	testIndices = indexer.GetSample(binSize, gen)
	testIndices = numpy.insert(testIndices, 0, [0, indexer.max() - 1])
	testIndices = numpy.unique(testIndices)
	
	sampleDeviations = list()
	binEdges = list()
	binWeights = list()
	for ii, i in enumerate(testIndices):
		for j in testIndices[ii + 1:]:
			if ((j - i) > 1):
				binEdges += [[sample[i], sample[j]],]
				binWeights += [integ.quad(dist.PDF, sample[i], sample[j])]
				sampleDeviations += [(binWeights[-1][0] * sampleSize/(j - i - 1) - 1.)*math.sqrt(j - i - 1),]
				if (binWeights[-1][1] == 0.):
					print(binEdges[-1])
					print(binWeights[-1])
	
	sampleDeviations = numpy.sort(sampleDeviations)
	sampleDeviations_2 = sampleDeviations**2
	
	report += Banner("Verify that the sample matches the PDF")
	report += "Bin a random sample with random bin edges, assume Poissonian counting errors,\n"
	report += "then compare empirical counts to the expected count from numerical integration of the PDF.\n\n"
	report += "    random sample size: {:}\n".format(int(sampleSize))
	report += "    number of random bins: {:}\n\n".format(len(sampleDeviations))
	report += "    random bin count deviations (in units of Poissonian error bars):\n"
	report += "        min ({0:.2f}), max ({1:.2f}), rms ({2:.2f})\n\n".format(
		sampleDeviations[0], sampleDeviations[-1], 
		math.sqrt(numpy.sum(sampleDeviations_2)/len(sampleDeviations_2)))
		
	meanError = mean if (dist.Mean() == 0.) else (mean/dist.Mean() - 1.)
	meanError *= math.sqrt(len(sample))
			
	report += "-"*80 + "\n\n"
	report += "Relative difference of sample mean/variance versus theory (in units of Poissonian error)\n\n"
	report += "    mean: {:.2f}    variance: {:.2f}\n\n".format(
		meanError, (variance/dist.Variance() - 1.)*math.sqrt(len(sample)))
		
	if hasattr(dist, 'CDF'):
		report += Banner("Verify that the CDF matches the PDF")
		
		report += "Plugging the smallest/largest variates into the CDF and CompCDF should return Order({:.3e})\n\n".format(2./sampleSize)
		report += "    CDF(smallest):    {:.1e},    1 - CDF(largest):      {:.1e}\n".format(dist.CDF(sample[0]), 1. - dist.CDF(sample[-1]))
		report += "    CompCDF(largest): {:.1e},    1 - CompCDF(smallest): {:.1e}\n".format(dist.CompCDF(sample[-1]), 1. - dist.CompCDF(sample[0]))
		
		PDF_CDF_deviations = list(map(lambda edges, weight:
			[(weight[0] - (dist.CDF(edges[1]) - dist.CDF(edges[0])))/weight[1],
			 (weight[0] - (dist.CompCDF(edges[0]) - dist.CompCDF(edges[1])))/weight[1]], 
			 binEdges, binWeights))
		
		PDF_CDF_deviations = numpy.transpose(numpy.asarray(PDF_CDF_deviations))
		PDF_CDF_deviations = numpy.asarray(list(map(numpy.sort, PDF_CDF_deviations)))
		
		PDF_CDF_deviations_2 = PDF_CDF_deviations**2
		
		report += "\n" + "-"*80 + "\n\n"
		report += "Using the same random bins as the PDF test, compare the integrated PDF\n"
		report += "to the difference in CDF across each bins.\n\n"
		report += "    CDF vs. PDF deviation (in units of the Order(epsilon) PDF integration error estimate):\n"
		report += "        min ({0:.2f}), max ({1:.2f}), rms ({2:.2f})\n\n".format(
			PDF_CDF_deviations[0][0], PDF_CDF_deviations[0][-1], 
			math.sqrt(numpy.sum(PDF_CDF_deviations_2[0])/len(PDF_CDF_deviations_2[0])))
		report += "    CompCDF vs. PDF deviation (in units of the Order(epsilon) PDF integration error estimate):\n"
		report += "        min ({0:.2f}), max ({1:.2f}), rms ({2:.2f})\n\n".format(
			PDF_CDF_deviations[1][0], PDF_CDF_deviations[1][-1], 
			math.sqrt(numpy.sum(PDF_CDF_deviations_2[1])/len(PDF_CDF_deviations_2[1])))
			
	if hasattr(dist, 'Q_small'):
		
		uSample = numpy.asarray([gen.HalfU_uneven() for __ in range(int(math.sqrt(sampleSize)))])
		Q_small_reverse = numpy.asarray(list(map(lambda u: dist.CDF(dist.Q_small(u)), uSample)))
		Q_large_reverse = numpy.asarray(list(map(lambda u: dist.CompCDF(dist.Q_large(u)), uSample)))
		
		Q_small_reverse = numpy.sort(Q_small_reverse/uSample - 1.)
		Q_large_reverse = numpy.sort(Q_large_reverse/uSample - 1.)
		
		Q_small_reverse_2 = Q_small_reverse**2
		Q_large_reverse_2 = Q_large_reverse**2
		
		report += Banner("Verify that the CDF reverse the quantile function (the two are self-consistent).")
		report += "    CDF reverses Q_small:\n"
		report += "        min ({0:.2e}), max ({1:.2e}), rms ({2:.2e})\n\n".format(
			Q_small_reverse[0], Q_small_reverse[-1],
			math.sqrt(numpy.sum(Q_small_reverse_2)/len(Q_small_reverse_2)))
		report += "    CompCDF reverses Q_large:\n"
		report += "        min ({0:.2e}), max ({1:.2e}), rms ({2:.2e})\n\n".format(
			Q_large_reverse[0], Q_large_reverse[-1],
			math.sqrt(numpy.sum(Q_large_reverse_2)/len(Q_large_reverse_2)))
		
	return report

assert(BinaryAccumulate(numpy.asarray(range(1000))) == (999.*1000.)/2.)

distVec = [pqr.uniform(-1., 100.),
	pqr.standard_normal(),
	pqr.normal(-1., math.pi),
	pqr.log_normal(1., 0.5),
	pqr.exponential(0.1),
	pqr.pareto(1., math.pi),
	pqr.weibull(0.1, 1.2),
	pqr.logistic(10., 3.),
	pqr.log_logistic(math.pi, 5.),
	pqr.gammaDist(10., math.pi)]
	
report = ""

for dist in distVec:
	print("Analyzing ... " + str(dist))
	report += TestDistribution(dist, 1e7, 1000)
	
with open("validation_report.txt", "w") as file:
	file.write(report)
