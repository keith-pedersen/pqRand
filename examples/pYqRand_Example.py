# Add the path that contains pYqRand.so, so we can import it
import sys
sys.path.append("../")

import pYqRand

# To draw random numbers, we need a pYqRand.engine to store the PRNG
# We pass a bool to the constructor of pYqRand.engine to inform it
# whether or not to do an initial, automatic seeding from 
# std::random_device (/dev/urandom on most GNU/LINUX).
# The default argument is True, so if we pass no arguments, 
# we default to doing the initial seed.

gen1 = pYqRand.engine() # equivalent to "gen1 = engine(True)"

# We can store gen1's seeded initial state to a file, for auditing/reuse.
gen1.WriteState("test.seed");

########################################################################

# We can now seed another gen2 from gen1's stored seed.
# We pass False during construction, to elude the initial seed,
# then manually call Seed_FromFile().
gen2 = pYqRand.engine(False) # If we forget to pass False, it only wastes time, it's not detrimental
gen2.Seed_FromFile("test.seed");

# We can test that the engines are equal by comparing their state-strings
print "\ngen2 == gen1, using gen1's seed file? \t\t" + str(gen2.GetState() == gen1.GetState())

########################################################################

# We can also seed another engine directly from a state-string, instead of a file
gen2state = gen2.GetState()

gen3 = pYqRand.engine(False)
gen3.Seed_FromString(gen2state); 

print "gen3 == gen2, using gen2's state string? \t" + str(gen3.GetState() == gen2.GetState())

########################################################################

# We can also directly copy the internal state of an engine
# This is much faster, because there is no parsing binary to ASCII and back
gen4 = pYqRand.engine(False)
gen4.Seed_FromEngine(gen3)

print "gen4 == gen3, using gen3 directly? \t\t" + str(gen4.GetState() == gen3.GetState()) + "\n"

########################################################################

# We can now test engines by plugging them through the same distribution, 
# and look for correlations between the random streams by plotting
# each stream on its own axis
import matplotlib.pyplot as plt

# Create a standard normal distribution (mu = 0, sigma = 1)
dist = pYqRand.standard_normal()

# We draw from a distribution by feeding it a pYqRand.engine as an argument
# This is the same API as C++ std::normal_distribution
print "One variate from a standard_normal: " + str(dist(gen4))

n = 10**4 # how large a sample size
pointArea = .01 # size of the plot point

# Right now, gen1 and gen2 are in the same state
# So their random streams should be totally correlated

fig1 = plt.figure(num="Use gen.Jump() to make independent equivalent generators independent")

fig1a = fig1.add_subplot(1,2,1)
fig1a.set_aspect('equal')
fig1a.set_xlabel('gen1')
fig1a.set_ylabel('gen2')
fig1a.set_title('Before Jump (correlated)')

x1 = [dist(gen1) for _ in range(n)]
y1 = [dist(gen2) for _ in range(n)]
fig1a.scatter(x1, y1, s = pointArea)

fig1b = fig1.add_subplot(1,2,2)
fig1b.set_aspect('equal')
fig1b.set_xlabel('gen1')
fig1b.set_ylabel('gen2')
fig1b.set_title('After Jump (uncorrelated)')

# Now we jump gen2, so it is no longer correlated to gen1
gen2.Jump()

x2 = [dist(gen1) for _ in range(n)]
y2 = [dist(gen2) for _ in range(n)]
fig1b.scatter(x2, y2, s = pointArea)

########################################################################

# Instead of having to jump manually, we can create a list of generators
# where the first is seeded from a given seed (using same conventions as Seed), 
# and the rest are jumped # once more than their left neighbor.
# We can test every pair of generators for correlations

fig2 = plt.figure(num="Use engine.ParallelList() to automatically make n independent generators")
fig2.subplots_adjust(wspace = .5, hspace = .5)
nGen = 4
subFigIndex = 1
engList = pYqRand.engine.ParallelList_FromFile(nGen, "test.seed")

# Record the initial seeds of each generator, so we can reset each to its
# initial state when we test each pair for correlations
engSeeds = [gen.GetState() for gen in engList]

for i in range(nGen):
	for j in range(i + 1, nGen):
		subfig = fig2.add_subplot(2,3,subFigIndex)
		subfig.set_aspect('equal')
		subfig.set_xlabel('gen' + str(i+1))
		subfig.set_ylabel('gen' + str(j+1))
		subFigIndex += 1
		gen_i = engList[i]
		gen_j = engList[j]
		# Reset each generator to its original state, so we're comparing
		# apples to apples
		gen_i.Seed_FromString(engSeeds[i])
		gen_j.Seed_FromString(engSeeds[j])
		x = [dist(gen_i) for _ in range(n)]
		y = [dist(gen_j) for _ in range(n)]
		subfig.scatter(x, y, s = pointArea)

########################################################################

# Bin the distribution into a histogram

nHist = 10**7
fdBins = int(14./(2.*0.937277/(float(nHist)**(1./3.)))) # Freedman-Draconis rule (not implemented in numpy 1.7.1)

fig3 = plt.figure(num="Bin standard_normal")
plt.hist([dist(gen3) for _ in range(nHist)], bins=fdBins, normed=1, histtype="step")

########################################################################

# Show all figures

plt.show()
