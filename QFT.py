#!/usr/bin/env python
# coding: utf-8

# In[14]:


import numpy as np
from numpy import pi
# importing Qiskit
from qiskit import QuantumCircuit, transpile, assemble, Aer, IBMQ
from qiskit.providers.ibmq import least_busy
from qiskit.tools.monitor import job_monitor
from qiskit.visualization import plot_histogram, plot_bloch_multivector

qc = QuantumCircuit(3)



# generalized version for N qbits
# this provides all the rotations for the circuit, does not do swaps
def qft_rotations(circuit, n):
    if n == 0: # if there are no qbits, dont do anything
        return circuit
    n -= 1 # n = n-1 because the circuit starts at 0 
    circuit.h(n) # hadamars the most significant qbit (the nth qbit)
    for qubit in range(n):
        # for each qbit above it, rotate it by the correct amount
        # the correct amoutn being pi/(2^x) where x is the degree of significance
        circuit.cp(pi/(2**(n-qubit)), qubit, n)
    # because we already reduced n by 1, to rotate all bits, run the program recusivly
    qft_rotations(circuit, n)

# propperly swaps registers
def swap_registers(circuit, n):
    for qubit in range (n//2):
        circuit.swap(qubit, n-1-qubit)
    return circuit



def qft(circuit, n):
    """QFT on the first N qubits in the circuit"""
    qft_rotations(circuit, n)
    swap_registers(circuit, n)
    return circuit



x = 4
cir = QuantumCircuit(x)
qft(cir,x)
cir.draw()



q = QuantumCircuit(3)
q.x(0)
q.x(2)
qft(q,3)
q.draw()



# make sure it does the thing
gem = QuantumCircuit(3)
# encode the binary number (5, 101)
# all qbits are 0 on creation so running them through an x gate switches them to 1
gem.x(0)
gem.x(2)
sim = Aer.get_backend("aer_simulator")
qc_init = gem.copy()
qc_init.save_statevector()
statevector = sim.run(qc_init).result().get_statevector()
plot_bloch_multivector(statevector)


# In[20]:


# run the qft on the circuit
qft(gem,3)
# create the thing
gem.save_statevector()
statevector = sim.run(gem).result().get_statevector()
plot_bloch_multivector(statevector)



def inverse_qft(circuit, n):
    # takes a number in fourier form and transforms it into classical bits
    # create a cirucit of the correct size and invert it
    invqft_circ = qft(QuantumCircuit(n),n).inverse()
    # stick it on the circuit containing our fourier number
    circuit.append(invqft_circ, circuit.qubits[:n])
    # decompose it and return
    # .decompose just lets us see the individual gates instead of a big U gate
    return circuit.decompose()


# create a circuit holding a fourier state n
nqubits = 3
number = 5
qc = QuantumCircuit(nqubits)
for qubit in range(nqubits):
    qc.h(qubit)
qc.p(number*pi/4,0)
qc.p(number*pi/2,1)
qc.p(number*pi,2)
# apply inverse qft
qc = inverse_qft(qc, nqubits)
qc.measure_all()
qc.draw()



IBMQ.load_account()
provider = IBMQ.get_provider(hub='ibm-q')
backend = least_busy(provider.backends(filters=lambda x: x.configuration().n_qubits >= nqubits 
                                       and not x.configuration().simulator 
                                       and x.status().operational==True))
print("least busy backend: ", backend)



shots = 2048
transpiled_qc = transpile(qc, backend, optimization_level=3)
job = backend.run(transpiled_qc, shots=shots)
job_monitor(job)


counts = job.result().get_counts()
plot_histogram(counts)

