import numpy as np
# Importing standard Qiskit libraries
from qiskit import QuantumCircuit, transpile, IBMQ
from qiskit.tools.jupyter import *
from qiskit.visualization import *
from ibm_quantum_widgets import *
from qiskit.visualization import plot_histogram, plot_bloch_multivector
from qiskit.providers.ibmq import least_busy
from qiskit.tools.monitor import job_monitor

Bits = 2
QC = QuantumCircuit(Bits,Bits)
QC.h(0)
QC.cnot(0,1)
QC.measure(range(Bits),range(Bits))
QC.draw()


IBMQ.load_account()
provider = IBMQ.get_provider(hub='ibm-q')
backend = least_busy(provider.backends(filters=lambda x: x.configuration().n_qubits >= Bits 
                                       and not x.configuration().simulator 
                                       and x.status().operational==True))
print("least busy backend: ", backend)


shots = 2048
transpiled_qc = transpile(QC, backend, optimization_level=3)
job = backend.run(transpiled_qc, shots=shots)
job_monitor(job)



counts = job.result().get_counts()
plot_histogram(counts)
