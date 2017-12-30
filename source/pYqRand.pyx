# This is a Cython wrapper for the pqRand C++ package
'''Randomly sample from important distributions (version 0.5.0)

	pYqRand generates random samples with the highest precision possible with floats

		import pYqRand as pqr

	The heart of the package is the engine storing the PRNG
	
		gen = pqr.engine()
	
	By default, an engine is automatically seeded from a good source of entropy, 
	but there is also a mechanism to store/retrieve (or manually specify) an engine's state.
	To draw a single variate, create a distribution object, then pass it an engine.
	
		expDist = pqr.exponential(1.)
		one = expDist(gen)
		many = expDist.GetSample(1e3, gen)

	pYqRand is a verbatim Cython wrapper of pqRand (C++), with a few minor differences:
	
	1. We cannot copy the state of a pYqRand.engine using assigment (operator=),
	because in Python this creates a reference, not a copy.
	To replace this functionality, we introduce ... pqr.engine.Seed_FromEngine
	
	2. We cannot pass a Python float by reference, so to give x a random sign:
		x = gen.ApplyRandomSign(x)
'''

########################################################################
########################################################################
########################################################################
# We first declare the C++ functions we intend to call
########################################################################

from libc.stdint cimport int64_t
from libcpp cimport bool # access bool as bool
from libcpp.string cimport string # access std::string as string
from libc.string cimport memcpy
from libcpp.vector cimport vector # access std::vector<T> as vector[T]
from cython.operator cimport dereference as deref 
# Dereference Cython pointers via deref() (e.g. func(deref(pointer)),
# to pass an object as argument. Note: pointer.func() has implicit dereference)

########################################################################

# According to Cython documentation, we need to manually replicate 
# typedefs declared in pqRand.hpp, but we merely need to get 
# the general type correct, without necessarily getting sizeof(type) exact.
# Thus, we will declare these typedefs using a type which is large enough 
# to accommodate either option in the C++ typedef.
ctypedef size_t result_type # could be uint32_t, but this is large enough
ctypedef double real_t # could be float/binary32, but this is large enough

########################################################################
# declare pqRand::engine
########################################################################

# If the pqRand::PRNG_t has a Jump() function, we can define engine's jump machinery.
# In spite of my best efforts, I cannot figure out how to set PRNG_CAN_JUMP based upon
# a compile-time constant in pqRand.hpp (e.g. its PRNG_CAN_JUMP).
# Hence, the following line must be set manually.
DEF PRNG_CAN_JUMP = 1

# It appears that you cannot add attributes to a cppclass after the original declaration.
# For example, the following will not work:
#
# 		cdef cppclass testClass
# 			bool func1()
#
# 		IF COMPILE_FLAG:
#			cdef cppclass testClass
#				bool func2()
#
# Thus, we need two declaration of engine_c, depending on PRNG_CAN_JUMP.
# This first (when PRNG_CAN_JUMP == True) will be fully commented.
# The second will merely remove the Jump() machinery from the first definition.
IF PRNG_CAN_JUMP:

	cdef extern from "pqRand/pqRand.hpp" namespace "pqRand":
		
		# To allow the objects to have the same name in Python and C++,
		# we declare the C++ classes as name_c, followed by the full C++ name in quotes.
		#	e.g. 		cdef nameInCython "nameInC++"
		cdef cppclass engine_c "pqRand::engine":
			
			# Functions declared as "except +" permit automatic conversion of their 
			# C++ exceptions into Python exceptions; otherwise,
			# thrown C++ exceptions cause a segmentation fault which crashes Python
			engine_c(const bool) except +
			
			# Cython doesn't support using "operator=" within the methods of the Python classes,
			# so we simply rename "operator=" for use in this file 
			engine_c& assign "pqRand::engine::operator=" (const engine_c&)
			
			result_type operator()()
					
			void Seed() except +
			# Note, cython does not recognize 'string const&', my preferred way to declare variables,
			# Similarly, "const" should be declared before a non-reference argument, 
			# otherwise I think that Cython is interpreting "const" as the name of the argument.
			void Seed_FromFile(const string&) except +
			void Seed_FromString(const string&) except +
			
			void WriteState(const string&) except +
			
			string GetState()
			
			bool RandBool()
			real_t ApplyRandomSign(const real_t)
			
			real_t U_uneven()
			real_t HalfU_uneven()
			real_t U_even()
			
			void Jump()
			vector[string] GetState_JumpVec(const size_t)
			
ELSE:
	cdef extern from "pqRand/pqRand.hpp" namespace "pqRand":
		
		cdef cppclass engine_c "pqRand::engine":
			
			engine_c(const bool) except +
			
			engine_c& assign "pqRand::engine::operator=" (const engine_c&)
			
			result_type operator()()
					
			void Seed() except +
			void Seed_FromFile(const string&) except +
			void Seed_FromString(const string&) except +
			
			void WriteState(const string&) except +
			
			string GetState()
			
			bool RandBool()
			real_t ApplyRandomSign(const real_t)
			
			real_t U_uneven()
			real_t HalfU_uneven()
			real_t U_even()
			
########################################################################
# declare distributions from distributions.hpp
########################################################################
		
cdef extern from "pqRand/distributions.hpp" namespace "pqRand":
	
	########################################################################
	# First declare the polymorphic interface classes
	########################################################################
	
	cdef cppclass distributionPDF_c "pqRand::distributionPDF":
		# min/max are not reserved words (like lambda)
		real_t min() const
		real_t max() const
		
		real_t operator()(engine_c& gen) const
		vector[real_t] GetSample(const size_t sampleSize, engine_c& gen) const
		
		real_t PDF(const real_t x) const
		real_t Mean() const
		real_t Variance() const
	
	cdef cppclass two:
		real_t x
		real_t y
		
	cdef two MeanAndVariance_c "pqRand::MeanAndVariance" (const distributionPDF_c& dist, 
		const size_t sampleSize, engine_c& gen)
		
	########################################################################
	
	cdef cppclass distributionCDF_c "pqRand::distributionCDF":
		real_t CDF(const real_t x) const
		real_t CompCDF(const real_t x) const
		
	########################################################################
	
	cdef cppclass distributionQ2_c "pqRand::distributionQ2":
		real_t Q_small(const real_t u) const
		real_t Q_large(const real_t u) const
			
	########################################################################
	# Now we declare the actual distributions, but only those functions
	# which are not already declared in the polymorphic parent classes
	########################################################################
	
	cdef cppclass uniform_c "pqRand::uniform":
		uniform_c(real_t const, real_t const) except +
	
	########################################################################
	
	cdef cppclass standard_normal_c "pqRand::standard_normal":
		standard_normal_c() except +
		
	########################################################################
	
	cdef cppclass normal_c "pqRand::normal":
		normal_c(real_t const, real_t const) except +

		real_t Mu() const
		real_t Sigma() const
		
	########################################################################
		
	cdef cppclass log_normal_c "pqRand::log_normal":
		log_normal_c(real_t const, real_t const) except +
		
		real_t Mu() const
		real_t Sigma() const	
		
	########################################################################
		
	cdef cppclass weibull_c "pqRand::weibull":
		weibull_c(real_t const, real_t const) except +

		real_t Lambda() const
		real_t k() const
		
	########################################################################
		
	cdef cppclass pareto_c "pqRand::pareto":
		pareto_c(real_t const, real_t const) except +

		real_t Alpha() const
		
	########################################################################
	
	cdef cppclass exponential_c "pqRand::exponential":
		exponential_c(real_t const) except +

		real_t Lambda() const
		
	########################################################################
	
	cdef cppclass logistic_c "pqRand::logistic":
		logistic_c(real_t const, real_t const) except +	
		
		real_t Mu() const
		real_t s() const
		
	########################################################################
	
	cdef cppclass log_logistic_c "pqRand::log_logistic":
		log_logistic_c(real_t const, real_t const) except +
		
		real_t Alpha() const
		real_t Beta() const
		
	########################################################################
	
	cdef cppclass gammaDist_c "pqRand::gammaDist":
		gammaDist_c(real_t const, real_t const) except +
		
		real_t Lambda() const
		real_t k() const
		
	########################################################################
	
	cdef cppclass uniform_integer_c "pqRand::uniform_integer<int64_t>":
		uniform_integer_c(const int64_t min_in, const int64_t max_in) except +
		
		int64_t operator()(engine_c& gen) const
		vector[int64_t] GetSample(const size_t sampleSize, engine_c& gen) const
		
		int64_t min() const
		int64_t max() const
		
########################################################################
########################################################################
########################################################################
# We now define the Cython wrapper classes of the C++ object.
########################################################################

########################################################################
# First we need a couple of helper functions, which are declared with a
# leading underscore to "hide" them from the user.

########################################################################		
# Python 3 uses unicode str; Python 2 uses ASCII strings; C++ std::string is ASCII only.
# Thus, Python 2 is a naive conversion between std::string and str, 
# whereas Python 3 needs to explicitly "encode"/"decode" the unicode.
# To allow this package to work with either Python version, we define two hidden helper functions
import sys

if(sys.version_info < (3,0)):
	def _str2string(s):
		return s
	def _string2str(s):
		return str(s)
else:
	def _str2string(s):
		return s.encode()
	def _string2str(s):
		return s.decode()
		
########################################################################		
# We cannot directly steal the data from a std::vector and 
# repackage it in a numpy.ndarray (like a std::move or swap),
# at least not without using undefined or implementation-specific behavior
# (hacking std::vector so it doesn't dellocate its data).
# Instead, we use a standard, quick copy. This is MUCH faster than 
# the automated version (return numpy.array(functionWhichReturnsStdVector()))

cimport numpy
import numpy

cdef _CopyToNumpy(const vector[real_t]& stdVec):
	# Using this "C API" for numpy is deprecated, as the compiler warns,
	# but I don't know any alternative
	cdef numpy.ndarray[double, ndim = 1, mode="c"] numpyVec = numpy.empty(stdVec.size())
	cdef double* vec = <double*>numpyVec.data
	
	# The memcpy appears to the be no faster than the loop, and is less readable
#~ 	memcpy(<void*>vec, <void*>stdVec.data(), stdVec.size()*sizeof(real_t))
	for i in range(stdVec.size()):
		vec[i] = stdVec[i]
	return numpyVec
	
cdef _CopyToNumpy_int(const vector[int64_t]& stdVec):
	# Using this "C API" for numpy is deprecated, as the compiler warns,
	# but I don't know any alternative
	cdef numpy.ndarray[int64_t, ndim = 1, mode="c"] numpyVec = numpy.empty(stdVec.size(), dtype=int) # This is the python int (aka int64_t)
	cdef int64_t* vec = <int64_t*>numpyVec.data
	
	# The memcpy appears to the be no faster than the loop, and is less readable
#~ 	memcpy(<void*>vec, <void*>stdVec.data(), stdVec.size()*sizeof(real_t))
	for i in range(stdVec.size()):
		vec[i] = stdVec[i]
	return numpyVec

########################################################################		
# First we declare the public pYqRand engine

# From <http://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html>
# "Cython initializes C++ class attributes of a cdef class using the nullary constructor. 
# If the class you’re wrapping does not have a nullary constructor, 
# you must store a pointer to the wrapped class and manually allocate and deallocate it. 
# A convenient and safe place to do so is in the __cinit__ and __dealloc__ methods 
# which are guaranteed to be called exactly once upon creation and deletion of the Python instance."

cdef class engine:
	'''
	gen = pqr.engine() ...  auto-seed the PRNG
	gen = pqr.engine(False) ... DO NOT seed the PRNG (e.g. to re-seed from stored).
	
	The engine object wraps the PRNG (xorshit1024*, period = 2**1024) with utility functions.
	The "Seed" functions allow storing states and re-seeding (or even manual seeding).
	In its simplest use, the user simply passes an engine to the various distributions,
	but they can also use the utility functions (e.g. uniform variates, random bool, random sign).'''
	cdef engine_c* c_engine
	
	def __cinit__(self, bool autoSeed = True):
		self.c_engine = new engine_c(autoSeed)
	
	def __dealloc__(self):
		del self.c_engine
	
	# operator()() in C++
	def __call__(self):
		'''
		gen() -> a random integer from [0,2**64)
		
		Calling the engine calls the underlying PRNG, returning a random, 64-bit, non-negative integer.'''
		return deref(self.c_engine)()
		
	def Seed(self):
		'''Auto-seed the PRNG from a good source of entropy (depending on the system, 
		either a hardware RNG or a cryptographic-PRNG seeded from environmental noise).'''
		self.c_engine.Seed()
		
	def Seed_FromFile(self, str filePath):
		'''Seed the PRNG from a file, interpreting the file's first line as a valid state-string
		
		Args:
			filePath (str): the location of the seed file
		
		Raises:
			OSError if the file cannot be found or opened.
			RuntimeError if the state-string is in the wrong format (see help(pqr.engine.GetState)).'''
		self.c_engine.Seed_FromFile(_str2string(filePath))
		
	def Seed_FromString(self, str state):
		'''
		Seed the PRNG from the supplied state (e.g. one obtained by calling gen.GetState()).
		
		Raises:
			RuntimeError if the state-string is in the wrong format
		
		Args:
			state (str): a specially formatted state-string (see engine.GetState())
			
		WARNING: 
			To manually seed the PRNG (i.e. from your own source of entropy), 
			ONLY SUPPLY the first 16 words of state (terminated with the number 16).
			The implementation will automatically (and repeatably) initialize the rest.
			Additionally, the PRNG has no "warm up", nor does it check that the seed is "good".
			If you supply a crappy seed (e.g. all zeroes), you will break PRNG uniformity.'''
		self.c_engine.Seed_FromString(_str2string(state))
	
	# Mimic C++ direct copy of engine state, a = b
	# The calling generator assumes the state of the argument generator.
	def Seed_FromEngine(self, engine original):
		'''Copy the internal state of the argument engine into this engine.'''
		self.c_engine.assign(deref(original.c_engine))
	
	def WriteState(self, str fileName):
		'''Write the PRNG's state-string to a file (overwriting it if it already exists).
						
		Args:
			filePath (str): the location of the seed file
		
		Raises:
			OSError if the file cannot be created or overwritten.'''
		self.c_engine.WriteState(_str2string(fileName))
		
	def GetState(self):
		'''Return the state of the PRNG as a state-string
		
		Format:			
			A valid state-string is a very specific single line containing
			human-readable, space-separated, non-negative integers:
			
			state_string = "word1␣word2␣...␣word16␣16␣p␣bitCache␣cacheMask"
			                                       ^^ this is the word count, so it is always 16
			
			 word[1-16] -- 16 space-separated, 64-bit integers from [0, 2**64)
			          p -- the index of the generator's active word, [0, 16)
			   bitCache -- a single 64-bit word whose bits become the random bools of RandBool()
			  cacheMask -- the location of the next bit used by RandBool
		'''
		return _string2str(self.c_engine.GetState())
		
	IF PRNG_CAN_JUMP:
		def Jump(self):
			'''Jump the generator forward by 2**512 calls (without actually calling it that many times).
			Useful for parallel threads, with independent engines, to prevent overlap in their sequences.'''
			self.c_engine.Jump()
			
		def GetState_JumpVec(self, int n):
			'''Given n, return a list of state-strings of length n, starting with the current state, 
			then calling Jump() once for each of the additional state.
			The PRNG ends in a state not in the list, so it is safe to keep using.'''
			return list(map(_string2str, self.c_engine.GetState_JumpVec(n)))
			
	def RandBool(self):
		'''Return a random bool (efficiently, using a cache of random bits, 1 bit per bool).'''
		return self.c_engine.RandBool()
	
	# Can't pass double by reference, so we return the altered victim
	def ApplyRandomSign(self, double victim):
		'''Return the argument (float) with a random sign.'''
		return self.c_engine.ApplyRandomSign(victim)
	
	def U_uneven(self):
		'''Return a random real number from the uniform distribution U(0,1], rounded to the nearest float.'''
		return self.c_engine.U_uneven()
		
	def HalfU_uneven(self): 
		'''Return a random real number from the uniform distribution U(0,0.5], rounded to the nearest float.'''
		return self.c_engine.HalfU_uneven()
		
	def U_even(self):
		'''Return a random real number from the uniform U(0,1), rounded to the nearest float in an
		EVENLY-spaced sample space (du = 2**-53).'''
		return self.c_engine.U_even()
		
########################################################################		
# Next we declare Python objects wrapping the polymorphic distribution classes
# These define the distribution interface, so that each Python object
# wrapping a given distribution does not need to do very much.

# We will keep the pointer to the C++ object in the most base object
cdef class _distributionPDF:
	cdef distributionPDF_c* dist
	
	# Construction and destruction must be handled by the final distribution class,
	# since all of the polymorphic classes are abstract.
	# We declare __cinit__ explicitly to avoid nullary ctor call from __init__
	def __cinit__(self):
		return
		
	def __dealloc__(self):
		del self.dist
		
	def min(self):
		'''The minimum variate sampled.'''
		return self.dist.min()

	def max(self):
		'''The maximum variate sampled.'''
		return self.dist.max()
	
	# IMPORTANT: even though this is a python function, using "float"
	# as x's type is not a Python float, but a C++ float (binary32).
	def PDF(self, double x):
		'''The probability distribution function (zero outside of [min, max]).'''
		return self.dist.PDF(x)
		
	def Mean(self):
		'''The distribution's mean.'''
		return self.dist.Mean()
	
	def Variance(self):
		'''The distribution's variance.'''
		return self.dist.Variance()
		
	def __call__(self, engine gen):
		'''
		Sample one variate by supplying a pYqRand.engine.
		
		Args:
			gen (engine): a pYqRand.engine PRNG'''
		return deref(self.dist)(deref(gen.c_engine))
	
	def GetSample(self, sampleSize, engine gen):
		'''
		Sample many variates and return them in a numpy.ndarray
		
		Args:
			  sampleSize: the number of variates to sample
			gen (engine): a pYqRand.engine PRNG
			
		Raises:
			ValueError if sampleSize is negative'''
		if (sampleSize < 0):
			raise ValueError("sampleSize must be non-negative")
		
		# type cast in Cython using carrots
		return _CopyToNumpy(self.dist.GetSample(<const size_t>sampleSize, deref(gen.c_engine)))
		
def MeanAndVariance(_distributionPDF distro, int sampleSize, engine gen):
	'''
	A validation function which calculates the mean and variance of a test sample, 
	for comparison to the analytic Mean() and Variance().
	
	Args:
		      distro: the distribution to test
		  sampleSize: the number of variates to sample
		gen (engine): a pYqRand.engine PRNG'''
	if(sampleSize < 0):
		raise ValueError("pYqRand::MeanAndVariance: sampleSize must be non-negative")
	
	meanVariance = MeanAndVariance_c(deref(distro.dist), <const size_t> sampleSize, deref(gen.c_engine))
	return [meanVariance.x, meanVariance.y]
		
########################################################################		

# To call the new functions available, we cast the base pointer
cdef class _distributionCDF(_distributionPDF):
	# Construction and destruction must be handled by the child class
	def __cinit__(self):
		return
		
	def CDF(self, double x):
		'''The cumulative distribution function.'''
		return  (<distributionCDF_c*> self.dist).CDF(x)
		
	def CompCDF(self, double x):
		'''The accurate complementary CDF (i.e. 1 - CDF(x), without cancellation).'''
		return  (<distributionCDF_c*> self.dist).CompCDF(x)
		
########################################################################
		
cdef class _distributionQ2(_distributionCDF):
	# Construction and destruction must be handled by the child class
	def __cinit__(self):
		return
		
	def Q_small(self, double u):
		'''The quantile function (the inverse of the CDF), which 
		accurately samples the small-value tail (given \p u < 1/2).'''
		return (<distributionQ2_c*> self.dist).Q_small(u)
		
	def Q_large(self, double u):
		'''The complementary quantile function 
		(the inverse of the Complementary CDF, by taking u -> 1 - u in Q_small),
		which accurately samples the large-value tail (given \p u < 1/2).'''
		return (<distributionQ2_c*> self.dist).Q_large(u)

########################################################################

cdef class uniform(_distributionCDF):
	'''
		uni = pqr.uniform(min, max)
	
	An object which samples from the uniform distribution [min, max], with PDF
	
		f(x) = (max-min)**-1      (min <= x <= max)
		
	Raises:
		ValueError if max <= min.'''	
		
	def __cinit__(self, double a, double b):
		self.dist = <distributionPDF_c*>(new uniform_c(a,b))
		
	def __str__(self):
		return "uniform distribution spanning the closed interval [{:.2e}, {:.2e}]".format(
			self.min(), self.max())

########################################################################

cdef class standard_normal(_distributionCDF):
	'''
		stdNorm = pqr.standard_normal()
	
	An object which samples from the standard normal distribution (mean = 0, variance = 1), with PDF
	
		f(x) = exp(-0.5 * x**2)/sqrt(2.*pi)'''
	
	def __cinit__(self):
		self.dist = <distributionPDF_c*>(new standard_normal_c())
		
	def __str__(self):
		return "standard normal distribution (a normal with mean Mu = 0 and std. dev. Sigma = 1)"
				
########################################################################
		
cdef class normal(_distributionCDF):
	'''
		norm = pqr.normal(Mu, Sigma)
	
	An object which samples from a normal distribution (Mu = mean, Sigma = std. dev.), with PDF
	
		f(x) = exp(-0.5 * (x-Mu)**2/Sigma**2)/sqrt(2.*pi*Sigma**2)
		
	Raises:
		ValueError is Sigma <= 0.'''
	
	def __cinit__(self, double mu, double sigma):
		self.dist = <distributionPDF_c*>(new normal_c(mu, sigma))
		
	def __str__(self):
		return "normal distribution with mean Mu = {:.2e} and std. dev. Sigma = {:.2e}".format(
			self.Mu(), self.Sigma())
		
	def Mu(self):
		'''The mean'''
		return (<normal_c*> self.dist).Mu()
		
	def Sigma(self):
		'''The standard deviation'''
		return (<normal_c*> self.dist).Sigma()
		
########################################################################
		
cdef class log_normal(_distributionCDF):
	'''
		logNorm = pqr.log_normal(Mu, Sigma)
	
	An object which samples from a log-normal distribution 
	(the log of variates have mean Mu and standard deviation Sigma), with PDF:
	
		f(x) = exp(-0.5 * (log(x)-Mu)**2/Sigma**2)/sqrt(2.*pi*Sigma**2)/x      (x >= 0)
		
	Raises:
		ValueError is Sigma <= 0.'''
	
	def __cinit__(self, double mu, double sigma):
		self.dist = <distributionPDF_c*>(new log_normal_c(mu, sigma))
		
	def __str__(self):
		return "log-normal distribution with Mu = {:.2e}, Sigma = {:.2e}".format(
			self.Mu(), self.Sigma())
		
	def Mu(self):
		'''The mean'''
		return (<log_normal_c*> self.dist).Mu()
		
	def Sigma(self):
		'''The standard deviation'''
		return (<log_normal_c*> self.dist).Sigma()		

########################################################################
		
cdef class weibull(_distributionQ2):
	'''
		weib = pqr.weibull(Lambda, k)
			
	An object which samples from a Weibull distribution (Lambda = scale, k = shape), with PDF
	
		f(x) = (k/Lambda) * (x/Lambda)**(k-1) * exp(-(x/Lambda)**k)      (x >= 0)
		
	Raises:
		ValueError if (Lambda <= 0) or (k <= 0).'''
	
	def __cinit__(self, double lambDUH, double k):
		self.dist = <distributionPDF_c*>(new weibull_c(lambDUH, k))
		
	def __str__(self):
		return "Weibull distribution with scale Lambda = {:.2e} and shape k = {:.2e}".format(
			self.Lambda(), self.k())

	def Lambda(self):
		'''The scale'''
		return (<weibull_c*> self.dist).Lambda()
		
	def k(self):
		'''The shape'''
		return (<weibull_c*> self.dist).k()

########################################################################
		
cdef class pareto(_distributionCDF):
	'''
		pareto = pqr.pareto(xMin, Alpha)
			
	An object which samples from a Pareto distribution (xMin = min, Alpha = index), with PDF
	
		f(x) = Alpha * xM**Alpha * x**-(Alpha + 1)      (x >= min)
		
	Raises:
		ValueError if (xMin <= 0) or (Alpha <= 0).'''
		
	def __cinit__(self, double xMin, double alpha):
		self.dist = <distributionPDF_c*>(new pareto_c(xMin, alpha))
	
	def __str__(self):
		return "Pareto distribution with minimum = {:.2e} and index Alpha = {:.2e}".format(
			self.min(), self.Alpha())

	def Alpha(self):
		'''The Pareto index'''
		return (<pareto_c*> self.dist).Alpha()

########################################################################

cdef class exponential(_distributionQ2):
	'''
		expo = pqr.exponential(Lambda)
			
	An object which samples from an exponential distribution (Lambda = rate), with PDF
	
		f(x) = Lambda * exp(-Lambda * x)      (x >= 0)
		
	Raises:
		ValueError if (Lambda <= 0).'''
	
	def __cinit__(self, double Lambda):
		self.dist = <distributionPDF_c*>(new exponential_c(Lambda))
		
	def __str__(self):
		return "exponential distribution with rate Lambda = {:.2e}".format(self.Lambda())
		
	def Lambda(self):
		'''The rate'''
		return (<exponential_c*> self.dist).Lambda()

########################################################################

cdef class logistic(_distributionQ2):
	'''
		logi = pqr.logistic(Mu, s)
			
	An object which samples from a logistic distribution with mean Mu and scale s, with PDF
	
		f(x) = cosh((x - Mu)/(2 s))**(-2) / (4 s)
		
	Raises:
		ValueError if (s <= 0).'''
	
	def __cinit__(self, double mu_in, double s_in):
		self.dist = <distributionPDF_c*>(new logistic_c(mu_in, s_in))
		
	def __str__(self):
		return "logistic distribution with mean Mu = {:.2e} and scale s = {:.2e}".format(
			self.Mu(), self.s())
		
	def Mu(self):
		'''The mean'''
		return (<logistic_c*> self.dist).Mu()
		
	def s(self):
		'''The scale'''
		return (<logistic_c*> self.dist).s()
		
########################################################################

cdef class log_logistic(_distributionQ2):
	'''
		ll = pqr.log_logistic(Alpha, Beta)
			
	An object which samples from a log-logistic distribution with scale Alpha and shape Beta, with PDF
	
		f(x) = (Beta/Alpha)*(x/Alpha)**(Beta-1) / (1 + (x/Alpha)**Beta)**2      (x >= 0)
		
	Raises:
		ValueError if (Alpha <= 0) or (Beta <= 0).'''
	
	def __cinit__(self, double alpha_in, double beta_in):
		self.dist = <distributionPDF_c*>(new log_logistic_c(alpha_in, beta_in))
		
	def __str__(self):
		return "log-logistic distribution with scale Alpha = {:.2e} and shape Beta = {:.2e}".format(
			self.Alpha(), self.Beta())
		
	def Alpha(self):
		'''The scale'''
		return (<log_logistic_c*> self.dist).Alpha()

	def Beta(self):
		'''The shape'''
		return (<log_logistic_c*> self.dist).Beta()

########################################################################	
		
cdef class gammaDist(_distributionPDF):
	'''
		gammaD = pqr.gammaDist(Lambda, k)
			
	An object which samples from a gamma distribution 
	(the sum of k exponential distributions with rate Lambda), with PDF
	
		f(x) = Lambda**k x**(k-1) exp(-Lambda * x)/Gamma(alpha)      (x >= 0)
		
	Rejection sampling is used, so it is not very fast.
		
	Raises:
		ValueError if (Lambda <= 0) or (k <= 0).'''
		
	def __cinit__(self, double lambda_in, double k_in):
		self.dist = <distributionPDF_c*>(new gammaDist_c(lambda_in, k_in))
		
	def __str__(self):
		return "gamma distribution with rate Lambda = {:.2e} and shape k = {:.2e}".format(
			self.Lambda(), self.k())
		
	def Lambda(self):
		'''The rate'''
		return (<gammaDist_c*> self.dist).Lambda()
		
	def k(self):
		'''The shape (how many exponentials are summed)'''
		return (<gammaDist_c*> self.dist).k()
		
cdef class uniform_integer:
	'''
		uni_int = pqr.uniform_integer(min, max)
			
	An object which samples integers from a uniform distribution 
	over the half-open interval [min, max)	(note that max is not sampled).
		
	Raises:
		ValueError if (max <= min).'''
	
	cdef uniform_integer_c* dist
	
	def __cinit__(self, int min_in, int max_in):
		self.dist = new uniform_integer_c(min_in, max_in)
		
	def __dealloc__(self):
		del self.dist
		
	def __str__(self):
		return "uniform integer distribution spanning the half-open interval [{}, {})".format(
			self.min(), self.max())
		
	def __call__(self, engine gen):
		'''
		Sample one variate by supplying a pYqRand.engine.
		
		Args:
			gen (engine): a pYqRand.engine PRNG'''
		return deref(self.dist)(deref(gen.c_engine))
	
	def GetSample(self, sampleSize, engine gen):
		'''
		Sample many variates and return them in a numpy.ndarray
		
		Args:
			  sampleSize: the number of variates to sample
			gen (engine): a pYqRand.engine PRNG
			
		Raises:
			ValueError if sampleSize is negative'''
		if (sampleSize < 0):
			raise ValueError("sampleSize must be non-negative")
		
		# type cast in Cython using carrots
		return _CopyToNumpy_int(self.dist.GetSample(<const size_t>sampleSize, deref(gen.c_engine)))
		
	def min(self):
		'''The minimum variate sampled.'''
		return self.dist.min()
		
	def max(self):
		'''One past the maximum variate sampled.'''
		return self.dist.max()
