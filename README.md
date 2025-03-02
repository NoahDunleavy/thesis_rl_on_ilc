The tracking of a system's output to some goal trajectory is a common industrial problem. Repetitive tasks conducted in a controlled manufacturing environment utilize complex machinery susceptible to noise and model characteristics not captured in their design or modelling process. Iterative Learning Control leverages this repetitious process in the presence of unknown, yet repeatable, disturbances to improve the output of each trial.

Constructing a controller which brings about this reduction in error is difficult to do without a system model. Reinforcement Learning helps overcome this by providing techniques to build controllers purely from input-output data. The number of data points needed to extract such a controller is the squared sum of the number of states and number of inputs.

When a system is translated into its ILC format, the number of effective states and inputs is scaled up by the number of steps in the manufacturing process. This exponentially increases the number of trials that would need to be run to produce a controller through RL. It is then desirable to reverse this increase in dimensions.

To accomplish this, we employ basis functions on the the system input and outputs. We find that the number of basis functions describing the input must be less than or equal to the number describing the output, and the input necessary to produce our desired output must be in the space of the input basis functions -- this requirement is not true for the output.

As we cannot know from the beginning what our goal input is, we must be able to dynamically grow our basis space representation in an efficient manner. To do so, we derive conjugate basis functions that are defined for a specific system.

We end with a methodology that allows for one to start with a low-dimension representation of a problem, learn the controller in the basis space through RL, and increase the dimensions as needed without compromising or changing the efficacy of previously learned parameters.
