# rice
A SPICE like simulator written entirely in Rust.

# Conventions

## Passive Sign Convention
For passive components (resistors, inductors, capacitors, etc...) positive current flows out of the positive voltage connection.

## Active Sign Convention
For active components (voltage sources, current sources, etc...) positive current flows into the positive voltage connection.

## Nodal Currents
For nodal currents positive current is current flowing out of the node and negative current is current flowing into the node.
