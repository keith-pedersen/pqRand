# distutils: language = c++
# distutils: sources = ['source/pqRand.cpp', 'source/distributions.cpp']
# distutils: extra_compile_args = ['--std=c++11']

# This is a Cython wrapper for the pqRand package
# It's annoying to have to write all these wrappers,
# but it's better than nothing

from libcpp cimport bool # access bool as bool
from libcpp.string cimport string # access std::string as string
# For fast cython construction, need pointer to ellude default initialization and seeding.
# Thus we need a way to dereference the pointer
from cython.operator cimport dereference as deref

########################################################################		

cdef extern from "include/pqRand.hpp" namespace "pqRand":
	cdef cppclass engine_c "pqRand::engine":
		engine_c(bool) except +
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
		
cdef extern from "include/distributions.hpp" namespace "pqRand":
	cdef cppclass standard_normal_c "pqRand::standard_normal":
		standard_normal_c() except +

		double operator()(engine_c& gen)
	
	cdef cppclass normal_c "pqRand::normal":
		normal_c(double const, double const) except +

		double operator()(engine_c& gen)
		
	cdef cppclass log_normal_c "pqRand::log_normal":
		log_normal_c(double const, double const) except +

		double operator()(engine_c& gen)
		
	cdef cppclass weibull_c "pqRand::weibull":
		weibull_c(double const, double const) except +

		double operator()(engine_c& gen)	
		
	cdef cppclass pareto_c "pqRand::pareto":
		pareto_c(double const, double const) except +

		double operator()(engine_c& gen)
	
	cdef cppclass exponential_c "pqRand::exponential":
		exponential_c(double const) except +

		double operator()(engine_c& gen)	

########################################################################		

cdef class engine:
	cdef engine_c* c_engine
	
	def __cinit__(self, bool doDefaultSeed = True):
		self.c_engine = new engine_c(doDefaultSeed)
	
	def __dealloc__(self):
		del self.c_engine
	
	def __call__(self):
		return deref(self.c_engine)()
		
	def Jump(self, unsigned int nTimes):
		deref(self.c_engine).Jump(nTimes)
		
	def Seed(self):
		deref(self.c_engine).Seed()
		
	def Seed_FromFile(self, str fileName):
		deref(self.c_engine).Seed_FromFile(fileName)
		
	def Seed_FromString(self, str seed):
		deref(self.c_engine).Seed_FromString(seed)
	
	def WriteState(self, str fileName):
		deref(self.c_engine).WriteState(fileName)
		
	def GetState(self):
		return deref(self.c_engine).GetState()
			
	def RandBool(self):
		return deref(self.c_engine).RandBool()
		
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
		
########################################################################
		
cdef class log_normal:
	cdef log_normal_c* dist
	
	def __cinit__(self, double mu, double sigma):
		self.dist = new log_normal_c(mu, sigma)
		
	def __dealloc__(self):
		del self.dist
	
	def __call__(self, engine gen):
		return deref(self.dist)(deref(gen.c_engine))

########################################################################
		
cdef class weibull:
	cdef weibull_c* dist
	
	def __cinit__(self, double lambDUH, double k):
		self.dist = new weibull_c(lambDUH, k)
		
	def __dealloc__(self):
		del self.dist
	
	def __call__(self, engine gen):
		return deref(self.dist)(deref(gen.c_engine))

########################################################################
		
cdef class pareto:
	cdef pareto_c* dist
	
	def __cinit__(self, double x_m, double alpha):
		self.dist = new pareto_c(x_m, alpha)
		
	def __dealloc__(self):
		del self.dist
	
	def __call__(self, engine gen):
		return deref(self.dist)(deref(gen.c_engine))
		
########################################################################
		
cdef class exponential:
	cdef exponential_c* dist
	
	def __cinit__(self, double mu):
		self.dist = new exponential_c(mu)
		
	def __dealloc__(self):
		del self.dist
	
	def __call__(self, engine gen):
		return deref(self.dist)(deref(gen.c_engine))
