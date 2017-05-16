# distutils: language = c++
# distutils: sources = ['source/pqRand.cpp', 'source/distributions.cpp']
# distutils: extra_compile_args = ['--std=c++11']

# This is a Cython wrapper for the pqRand C++ package
# It's annoying to have to write all these wrappers,
# but it's better than nothing.

# Differences versus C++ package
# 	1. We cannot copy the state of an engine using assigment (operator=),
# 		because in Python this creates a reference, not a copy.
# 		To replace this functionality, we introduce the function Seed_FromEngine().
#  2. engine.ApplyRandomSign(val) cannot alter val, so we return the altered val instead.
#  3. Added ParallelList functions --- test to see if this is useful in Python before moving to C++

from libcpp cimport bool # access bool as bool
from libcpp.string cimport string # access std::string as string
# For fast cython construction, need pointer to ellude default initialization and seeding.
# Thus we need a way to dereference the pointer
from cython.operator cimport dereference as deref

########################################################################
# First we must declare the c-functions we intend to call

########################################################################
# import pqRand.hpp

# To allow the objects to have the same name in Cython and C++,
# we import as **_c, giving the full C++ name in quotes
cdef extern from "include/pqRand.hpp" namespace "pqRand":
	cdef cppclass engine_c "pqRand::engine":
		engine_c(bool) except +
		
		# Cython doesn't support access to operator=, 
		# so we simply rename the functionality
		engine_c& assign "pqRand::engine::operator=" (engine_c& const)
		
		unsigned long operator()()
		void Jump(unsigned long nTimes)
		
		void Seed()
		# Note, cython does not recognize 'string const&', my preferred way to declare this variable 
		void Seed_FromFile(string& const)
		void Seed_FromString(string& const)
		
		void WriteState(string& const)
		string GetState()
		
		bool RandBool()
		void ApplyRandomSign(double&)
		
		double U_Q()
		double HalfU_Q()
		double U_S()
		double HalfU_S()
		double U_S_canonical()
		
########################################################################		
# import distributions.hpp
		
cdef extern from "include/distributions.hpp" namespace "pqRand":
	cdef cppclass standard_normal_c "pqRand::standard_normal":
		standard_normal_c() except +

		double operator()(engine_c& gen)
	
	cdef cppclass normal_c "pqRand::normal":
		normal_c(double const, double const) except +

		double operator()(engine_c& gen)
		double Mu()
		double Sigma()
		
	cdef cppclass log_normal_c "pqRand::log_normal":
		log_normal_c(double const, double const) except +
		double Mu()
		double Sigma()

		double operator()(engine_c& gen)
		
	cdef cppclass weibull_c "pqRand::weibull":
		weibull_c(double const, double const) except +

		double operator()(engine_c& gen)
		double Lambda()
		double k()
		
	cdef cppclass pareto_c "pqRand::pareto":
		pareto_c(double const, double const) except +

		double operator()(engine_c& gen)
		double xM()
		double Alpha()
	
	cdef cppclass exponential_c "pqRand::exponential":
		exponential_c(double const) except +

		double operator()(engine_c& gen)
		double Lambda()

########################################################################		
# Now we make the Cython wrapper class of the C++ object

# Cython default initializes the C++ object using the default ctor
# Howerver, engine() causes the default seeding, which we do not always want.
# To ellude the auto-default behavior, we must define __cinit__ and __dealloc__, 
# which in this case requires storing a pointer to the C++ object

cdef class engine:
	cdef engine_c* c_engine
	
	def __cinit__(self, bool doDefaultSeed = True):
		self.c_engine = new engine_c(doDefaultSeed)
	
	def __dealloc__(self):
		del self.c_engine
	
	# operator()() in C++
	def __call__(self):
		return deref(self.c_engine)()
	
	# No overloading in Python, use default argument to mimic C++ behavior
	def Jump(self, unsigned int nTimes = 1):
		deref(self.c_engine).Jump(nTimes)
	
	def Seed(self):
		deref(self.c_engine).Seed()
		
	def Seed_FromFile(self, str fileName):
		deref(self.c_engine).Seed_FromFile(fileName)
		
	def Seed_FromString(self, str seed):
		deref(self.c_engine).Seed_FromString(seed)
	
	# Mimic C++ direct copy of engine state, a = b
	# The calling generator assumes the state of the argument generator.
	def Seed_FromEngine(self, engine original):
		deref(self.c_engine).assign(deref(original.c_engine))
	
	def WriteState(self, str fileName):
		deref(self.c_engine).WriteState(fileName)
		
	def GetState(self):
		return deref(self.c_engine).GetState()
			
	def RandBool(self):
		return deref(self.c_engine).RandBool()
	
	# Can't pass double by reference, so we return the altered victim
	def ApplyRandomSign(self, double victim):
		deref(self.c_engine).ApplyRandomSign(victim)
		return victim
	
	def U_Q(self): 
		return deref(self.c_engine).U_Q()
		
	def HalfU_Q(self): 
		return deref(self.c_engine).HalfU_Q()
		
	def U_S(self): 
		return deref(self.c_engine).U_S()
		
	def HalfU_S(self): 
		return deref(self.c_engine).HalfU_S()
		
	def U_S_canonical(self):
		return deref(self.c_engine).U_S_canonical()
	
	# Introduce helper funtions to create parallel states of the generator.
	
	# The worker function used by the three "public" functions.
	# Generate nEngines in orthogonal state by seeding the first generator,
	# then jumping the rest of the generators.	
	# Make the first index generator equal to the passed generator, 
	# then Jump() once for each additional additional index
	@staticmethod
	def __ParallelList_FromEngine(unsigned int nEngines, engine first):
		engList = [engine(False) for _ in range(nEngines)]
		engList[0].Seed_FromEngine(first)
		for i in range(1, nEngines):
			engList[i].Seed_FromEngine(engList[i-1])
			engList[i].Jump()
		return engList
	
	@staticmethod
	def ParallelList(unsigned int nEngines):
		first = engine() # auto-seed
		return engine.__ParallelList_FromEngine(nEngines, first)
	
	@staticmethod
	def ParallelList_FromFile(unsigned int nEngines, str fileName):
		first = engine(False)
		first.Seed_FromFile(fileName)
		return engine.__ParallelList_FromEngine(nEngines, first)
		
	@staticmethod
	def ParallelList_FromString(unsigned int nEngines, str seed):
		first = engine(False)
		first.Seed_FromString(seed)
		return engine.__ParallelList_FromEngine(nEngines, first)	
		
########################################################################

cdef class standard_normal:
	cdef standard_normal_c* dist
	
	def __cinit__(self):
		self.dist = new standard_normal_c()
		
	def __dealloc__(self):
		del self.dist
	
	def __call__(self, engine gen):
		return deref(self.dist)(deref(gen.c_engine))
		
########################################################################
		
cdef class normal:
	cdef normal_c* dist
	
	def __cinit__(self, double mu, double sigma):
		self.dist = new normal_c(mu, sigma)
		
	def __dealloc__(self):
		del self.dist
	
	def __call__(self, engine gen):
		return deref(self.dist)(deref(gen.c_engine))
		
	def Mu(self):
		return deref(self.dist).Mu()
		
	def Sigma(self):
		return deref(self.dist).Sigma()
		
########################################################################
		
cdef class log_normal:
	cdef log_normal_c* dist
	
	def __cinit__(self, double mu, double sigma):
		self.dist = new log_normal_c(mu, sigma)
		
	def __dealloc__(self):
		del self.dist
	
	def __call__(self, engine gen):
		return deref(self.dist)(deref(gen.c_engine))
		
	def Mu(self):
		return deref(self.dist).Mu()
		
	def Sigma(self):
		return deref(self.dist).Sigma()		

########################################################################
		
cdef class weibull:
	cdef weibull_c* dist
	
	def __cinit__(self, double lambDUH, double k):
		self.dist = new weibull_c(lambDUH, k)
		
	def __dealloc__(self):
		del self.dist
	
	def __call__(self, engine gen):
		return deref(self.dist)(deref(gen.c_engine))

	def Lambda(self):
		return deref(self.dist).Lambda()
		
	def k(self):
		return deref(self.dist).k()

########################################################################
		
cdef class pareto:
	cdef pareto_c* dist
	
	def __cinit__(self, double x_m, double alpha):
		self.dist = new pareto_c(x_m, alpha)
		
	def __dealloc__(self):
		del self.dist

	def __call__(self, engine gen):
		return deref(self.dist)(deref(gen.c_engine))

	def xM(self):
		return deref(self.dist).xM()

	def Alpha(self):
		return deref(self.dist).Alpha()

########################################################################
		
cdef class exponential:
	cdef exponential_c* dist
	
	def __cinit__(self, double Lambda):
		self.dist = new exponential_c(Lambda)
		
	def __dealloc__(self):
		del self.dist
	
	def __call__(self, engine gen):
		return deref(self.dist)(deref(gen.c_engine))
		
	def Lambda(self):
		return deref(self.dist).Lambda()
