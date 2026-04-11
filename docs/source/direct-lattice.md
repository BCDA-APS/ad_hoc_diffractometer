# Direct lattice in crystallography

## Vector math order of operations

> Given the vector math expression: A dot B cross C, what is the order of operations?

In the vector math expression **A · (B × C)**, the **cross product has
higher precedence than the dot product**, so the order of operations is:

1. **First: Calculate B × C** (the cross product)
2. **Second: Calculate A · (result)** (the dot product)

This is the standard convention in vector mathematics. The cross product
is always evaluated before the dot product unless parentheses explicitly
indicate otherwise.

### Why this order matters

This particular combination—**A · (B × C)**—is called the **scalar
triple product** (or mixed product). It has a special geometric meaning:
it gives the signed volume of the parallelepiped formed by the three
vectors A, B, and C.

If you wanted to reverse the order, you'd need to write it explicitly as
**(A · B) × C**, which would give you a vector result instead of a
scalar, and it would have a completely different meaning.

------------------------

I need to clarify the representation you're asking about, because there
are different conventions for encoding crystal lattice parameters into
vectors.

The most common approach in crystallography is to represent the lattice
parameters using **lattice vectors** (the edge vectors of the unit
cell). In this case:

| Vector | Represents | Formula |
|--------|-----------|---------|
| **A** | Lattice vector a | **(a, 0, 0)** |
| **B** | Lattice vector b | **(b·cos(γ), b·sin(γ), 0)** |
| **C** | Lattice vector c | **(c·cos(β), c·sin(α)·sin(β), c·cos(α))** |

where:

- **a, b, c** are the lattice parameters (edge lengths)
- **α, β, γ** are the angles between the edges (α between b and c, β
  between a and c, γ between a and b)

### Why this representation?

With vectors defined this way, **A · (B × C)** gives you the **unit cell
volume**—a fundamental property in crystallography. This is exactly why
the scalar triple product is so important in crystallography
calculations.

## Standard form of lattice

Perfect. Here are the **lattice vectors** in standard crystallographic form:

| Vector | Formula |
|--------|---------|
| **A** | **(a, 0, 0)** |
| **B** | **(b·cos(γ), b·sin(γ), 0)** |
| **C** | **(c·cos(β), c·sin(α)·sin(β), c·cos(α))** |

### Breaking down each vector

**Vector A** is the simplest—it lies along the x-axis:

- Component 1: **a**
- Component 2: **0**
- Component 3: **0**

**Vector B** lies in the xy-plane, making angle **γ** (gamma) with vector A:

- Component 1: **b·cos(γ)**
- Component 2: **b·sin(γ)**
- Component 3: **0**

**Vector C** points into 3D space, with **β** (beta) defining its angle
from the x-axis and **α** (alpha) defining the tilt:

- Component 1: **c·cos(β)**
- Component 2: **c·sin(α)·sin(β)**
- Component 3: **c·cos(α)**

### Key point

When you compute **A · (B × C)** using these vectors, you get the **unit
cell volume** of the crystal lattice. This is the fundamental reason
this scalar triple product is so important in crystallography.
