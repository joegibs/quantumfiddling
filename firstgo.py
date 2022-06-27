# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 17:18:39 2022

@author: jogib
"""

import cirq
import numpy as np

#%%
qubit = cirq.GridQubit(0, 0)

circuit = cirq.Circuit(cirq.X(qubit) ** 0.5, cirq.measure(qubit, key="m"))


print("Circuit:\n", circuit)

simulator = cirq.Simulator()
result = simulator.run(circuit, repetitions=30)

print("Results:\n", result)

#%%

a = cirq.NamedQubit("a")

ops = [cirq.rx(np.pi / 2).on(a), cirq.measure(a, key="m")]
circuit = cirq.Circuit(ops)
print("Circuit:\n", circuit)
simulator = cirq.Simulator()
result = simulator.run(circuit, repetitions=30000)
# _ = cirq.vis.plot_state_histogram(result, plt.subplot())
# print('Results:\n',result)

state_vector_before_measurement = simulator.simulate(circuit[:-1])
sampled_state_vector_after_measurement = simulator.simulate(circuit)

print(f"State before measurement:")
print(state_vector_before_measurement)
print(f"State after measurement:")
print(sampled_state_vector_after_measurement)
