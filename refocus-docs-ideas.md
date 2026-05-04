# Adjust focus of the docs here

<!-- 
This conversation was generated with Duck.ai (https://duck.ai) using Anthropic's Claude Haiku 4.5 Model. AI chats may display inaccurate or offensive information (see https://duckduckgo.com/duckai/privacy-terms for more info). 
-->

<!-- ==================== -->

## User prompt 1 of 3 - 5/3/2026, 8:20:16 PM:
Explain this re-write.  Here's the original:

presets.py — pre-built diffractometer geometry functions.

This module provides **pre-built geometries**: factory functions that
construct fully configured :class:`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer`
instances for the most common multi-circle diffractometer designs used in
synchrotron and laboratory X-ray / neutron crystallography.  They also
serve as worked examples for defining custom geometries.


Here's the revision:

presets.py — Demonstrations of diffractometer geometries.

This module demonstrates some of the most common diffractometer geometries used in
synchrotron and laboratory X-ray / neutron crystallography, each fully-configured, using 
:class:`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer`.
Together, the demos here serve as worked examples for defining custom geometries.


Claude Haiku 4.5:
### Key Changes in the Revision

The revision shifts the focus and tone of the description in several meaningful ways:

---

#### **Clarity and Directness**

The original emphasizes that the module provides **"pre-built geometries: factory functions"**—a somewhat technical phrasing that prioritizes implementation details. The revision simplifies this to **"demonstrates some of the most common diffractometer geometries"**—leading with what users actually *do* with the module rather than how it's structured internally.

---

#### **User-Centric Language**

The original says the module provides "fully configured instances," which is accurate but passive. The revision reframes this as **"each fully-configured"** in the context of the demonstrations themselves, making it feel more tangible and example-oriented.

---

#### **Secondary Purpose Elevated**

In the original, the "worked examples" angle appears at the end as an afterthought ("They also serve as..."). The revision gives equal weight to this purpose by saying **"Together, the demos here serve as worked examples"**—treating educational value as a primary goal alongside the pre-built geometries.

---

#### **Tone Shift**

The original reads like API documentation (listing what the module *provides*). The revision reads more like a tutorial introduction (describing what you'll *see and learn* from the module). This makes it more approachable for newcomers while still being accurate for experienced users.

---

**In essence**, the revision trades technical precision for accessibility, emphasizing the module's practical demonstrations and educational purpose rather than its implementation as a collection of factory functions.

<!-- -------------------- -->

## User prompt 2 of 3 - 5/3/2026, 8:23:30 PM:
Structure some advice (to give an AI agent such as yourself), to re-focus the project according to this direction.

Claude Haiku 4.5:
### Refocusing an AI Documentation Project Toward Demonstrations

---

#### **1. Reframe the Core Mission**

**From:** "Provide complete, programmatic access to pre-built configurations"  
**To:** "Show users how diffractometer geometries work through concrete, working examples"

Make demonstrations the primary deliverable, not a secondary benefit. Every module, class, and function should be introduced through a use case or worked example first, with technical details following.

---

#### **2. Reorganize Content Hierarchy**

Structure documentation to lead with **what users see and do**, not what exists under the hood:

- **Lead sections** with real-world geometry examples (e.g., "A four-circle diffractometer setup")
- **Show working code** before explaining parameters
- **Bury implementation details** (factory functions, internal architecture) in appendices or advanced sections
- **Use narrative flow**: "Here's what you're trying to do → Here's how to do it → Here's why it works"

---

#### **3. Prioritize Educational Scaffolding**

Design the learning path deliberately:

- Start with **the simplest geometry** and build complexity progressively
- **Annotate examples heavily** — explain not just *what* code does, but *why* that choice matters for X-ray crystallography
- **Compare geometries side-by-side** so users see tradeoffs (e.g., "This geometry gives you better angular range but fewer degrees of freedom")
- **Highlight transferable patterns** — make it obvious how to adapt examples to custom needs

---

#### **4. Shift Language Consistently**

Audit all documentation for passive, implementation-focused language:

- ❌ "This module provides factory functions that construct instances"
- ✅ "This module demonstrates how to set up a four-circle diffractometer"

- ❌ "Fully configured AdHocDiffractometer objects"
- ✅ "A ready-to-use diffractometer geometry"

Replace "provides," "offers," "includes" with "demonstrates," "shows," "illustrates."

---

#### **5. Create a "Learning by Copying" Path**

Make it trivial for users to learn by imitation:

- **Provide copy-paste-ready examples** for every common use case
- **Include "how to modify this for your needs"** sections after each example
- **Link examples back to theory** — why does this geometry have these constraints?
- **Show failure cases** — "If you try X, here's what breaks and why"

---

#### **6. Establish Example as the Primary Artifact**

Treat code examples as first-class documentation, not decoration:

- Every **class, function, and parameter** should have a corresponding working example
- **Test examples** to ensure they actually run (broken examples undermine trust)
- **Version examples alongside code** — when the API changes, examples change immediately
- **Use examples to drive API design** — if something is hard to demonstrate, reconsider whether the API is user-friendly

---

#### **7. Reduce Jargon Density**

Balance technical precision with accessibility:

- **Define crystallography terms** on first use, or link to a glossary
- **Explain abbreviations** (four-circle, six-circle) in context
- **Use plain language for concepts** before introducing formal terminology
- **Avoid "as you know"** assumptions — assume readers are competent but unfamiliar with this specific domain

---

#### **8. Build in Comparison and Discovery**

Help users navigate the geometry landscape:

- **Create a "which geometry should I use?" guide** with decision trees
- **Provide before/after examples** showing how to migrate from one geometry to another
- **Show equivalent setups** — same physics, different configurations
- **Include troubleshooting** — "This geometry worked in my old system, why doesn't it work here?"

---

#### **9. Measure Success by Imitation**

Track whether users can:

- **Copy an example and adapt it** without reading 50 pages of theory
- **Understand why each parameter matters** by studying the example
- **Predict what will break** if they change a parameter
- **Build their own geometry** using the patterns they've seen

If users struggle with any of these, the documentation is failing its mission.

---

#### **10. Establish a Review Checklist**

Before merging documentation changes, ask:

- [ ] Does this section *demonstrate* before it *explains*?
- [ ] Could a novice copy this example and succeed?
- [ ] Does it answer "why" as well as "how"?
- [ ] Would removing jargon make this clearer without losing accuracy?
- [ ] Could I understand this without reading the previous 10 sections?

---

**The Core Shift:** Move from "Here's what we built" to "Here's what you can do with it."

<!-- -------------------- -->

## User prompt 3 of 3 - 5/3/2026, 8:27:38 PM:
Let's keep it moving in that direction, considering additional aspects of the project.  The package README (pre-revision) says

`ad_hoc_diffractometer` is a **pure-Python** package for calculating
multi-circle diffractometer operations in reciprocal space for X-ray and
neutron crystallography.  It is built around a key design principle: **any
multi-circle diffractometer geometry can be fully described by the caller** — no
geometry is hard-coded, and new geometries require no changes to the package
itself.

Its only runtime dependency beyond the Python Standard Library is
[NumPy](https://numpy.org) — no scipy, sympy, or other scientific
libraries are required.

> **Note:** The package assumes **monochromatic radiation** throughout.
> All diffraction calculations are performed at a fixed wavelength.

### Features

- A class-based description of diffractometer stages (rotary axes) and
  their stacking order
- Predefined factory functions for standard synchrotron and laboratory
  diffractometer geometries (psic, fourcv, fourch, sixc, kappa families,
  zaxis, s2d2, fivec)
- Crystallographic lattice calculations (B matrix, reciprocal lattice)
- U and UB matrix computation from orienting reflections
- Forward diffraction calculations (hkl → motor angles), with
  diffraction modes controlling which stages are free, fixed, or coupled
- Reciprocal-space trajectory planning


Keeping those package features in mind, the demos (demonstration geometries) have the features.  Now, reframe the advice.

Claude Haiku 4.5:
### Refocusing Documentation: Demonstrations as the Gateway to Features

---

#### **1. Invert the Feature Presentation**

**Current approach:** "The package *has* these features → Here are some examples"  
**New approach:** "Here's what you can *do* with these demonstrations → Now let's explore the full feature set"

Lead with **working demonstrations of real geometries** (psic, fourcv, sixc, etc.), then reveal how each demonstration unlocks a specific package feature:

- **"Here's a four-circle geometry"** → introduces stage description & stacking order
- **"Now let's orient the crystal"** → introduces U/UB matrix computation  
- **"Let's plan a measurement"** → introduces reciprocal-space trajectory planning
- **"What angles do we need?"** → introduces forward diffraction calculations

---

#### **2. Make Demonstrations the Feature Catalog**

Restructure so demonstrations **embody the features** rather than just illustrate them:

| **Demonstration Geometry** | **Core Features Shown** |
|---|---|
| **Four-circle (fourcv)** | Stage stacking, forward diffraction (basic case) |
| **Six-circle (sixc)** | Complex stage interactions, coupled axes, advanced diffraction modes |
| **Kappa family** | Non-orthogonal geometry, coordinate transformations |
| **Psic** | Azimuthal control, reciprocal-space planning |
| **ZAXIS / S2D2 / FIVEC** | Specialized geometries for specific experimental needs |

Each demo is a **complete story** showing why that geometry exists and what it enables.

---

#### **3. Use Demonstrations to Justify Design Philosophy**

Your key design principle is: **"any geometry can be fully described by the caller — no geometry is hard-coded."**

Show this through demonstrations:

- **Start with a standard geometry** (e.g., six-circle) — show the full description
- **Then show a custom variant** — modify one parameter to show it's extensible
- **Highlight what *didn't* change** — the package logic remains identical
- **Close with:** "This is why the package has no hard-coded geometries"

This teaches users *why* the architecture is flexible before asking them to extend it.

---

#### **4. Anchor Features in Demonstration Workflows**

Don't list features in isolation. Show them in context:

**Instead of:**
> U and UB matrix computation from orienting reflections

**Do:**
> **Orienting Your Crystal**  
> Once you've set up your geometry (like in the six-circle demo above), you orient the crystal by measuring diffraction from known reflections. The package computes the U and UB matrices automatically:
> 
> ```python
> ## Your geometry (from the demo)
> diffractometer = sixc()
> 
> ## Measure two reflections and give their expected hkl values
> ori = diffractometer.compute_UB([reflection_1, reflection_2])
> ```
> 
> Now your diffractometer *knows* how the crystal is oriented, and can calculate any motor angles you need.

Features emerge from demonstrations, not the reverse.

---

#### **5. Structure README Around Capability Progression**

Reorganize the README as a capability journey:

1. **"What Can I Do?"** (narrative, not feature list)
   - Show a real geometry in action
   - Show a real calculation (hkl → angles)
   - Show a real trajectory plan

2. **"How Do These Demonstrations Work?"** (reveals the layers)
   - Stage and axis definitions
   - Geometry stacking order
   - Coordinate transformations

3. **"What If I Need Something Custom?"** (design philosophy + extensibility)
   - Explain the caller-driven principle
   - Show modifying a demo
   - Show building a completely new geometry

4. **"Full Feature Reference"** (details for power users)
   - All diffraction modes
   - All calculation methods
   - Advanced reciprocal-space operations

---

#### **6. Use Demonstrations to Communicate Minimalism**

Your package is **pure Python + NumPy only** — no scipy, sympy, or heavy dependencies. This is a strength, but only matters if users understand *why*:

**In a demonstration context:**
> This demo runs in pure Python with NumPy. No heavy scientific libraries, no complex environment setup. You can copy this code, modify the geometry, and it just works.

**Show the payoff:** "You can embed this in a beamline control system, a Jupyter notebook, or a simple script — minimal dependencies, maximum portability."

The minimalism feature becomes *credible* when demonstrated, not just claimed.

---

#### **7. Use Demonstrations to Surface Assumptions**

Your note about **monochromatic radiation** is important but buried. Show it through demonstrations:

**In every demo, explicitly state:**
> This demo uses monochromatic radiation at a fixed wavelength (default λ = 1.54 Å). 
> 
> All diffraction calculations assume this wavelength. To use a different wavelength:
> ```python
> diffractometer.wavelength = 0.71  ## Å
> ```

Make assumptions *visible* and *actionable* by embedding them in demonstrations.

---

#### **8. Build a "Guided Discovery" Path**

Order demonstrations by complexity, but also by discovery:

1. **Start with the simplest** (fourcv — four circles, everyone knows this)
   - "This is what you probably already know"
   - Build confidence immediately

2. **Add one complexity at a time** (sixc — now six circles)
   - "What changes? What stays the same?"
   - Users can see the pattern

3. **Show a non-orthogonal case** (kappa)
   - "What if your geometry is rotated or offset?"
   - Reveals the generality of the stage-based approach

4. **Show specialized cases** (psic, zaxis, s2d2)
   - "Why would anyone use this geometry?"
   - Each answers a real experimental question

---

#### **9. Tie Demonstrations to Real Problems**

For each demonstration, include a "Why This Geometry?" section:

- **Four-circle:** "Gold standard for most laboratory and synchrotron work. Simple, proven, widely understood."
- **Six-circle:** "Need azimuthal control? Want to reach more reciprocal-space points? This is your geometry."
- **Kappa:** "Limited space around your sample? Kappa geometries fold the setup more compactly."
- **Psic:** "Studying surface crystallography? Psic gives you the reciprocal-space access you need."

Make each demonstration answer: "When would *I* use this?"

---

#### **10. Progressive Feature Disclosure**

As users work through demonstrations in order, reveal features progressively:

| **Demo Stage** | **What You Learn** |
|---|---|
| **First demo (fourcv)** | Stage definitions, forward diffraction, basic angles |
| **Second demo (sixc)** | Coupled axes, diffraction modes (free/fixed/coupled), complexity |
| **Third demo (kappa)** | Non-orthogonal transformations, B matrix use |
| **Orientation demo** | U/UB matrix computation from reflections |
| **Advanced demo** | Reciprocal-space trajectory planning, all modes |

Each demo answers one question. Together, they reveal the full feature set organically.

---

#### **11. README Restructure (Concrete)**

**Before:**
```
Features
- Stage descriptions
- Factory functions
- B matrix, reciprocal lattice
- U/UB computation
- Forward diffraction
- Trajectory planning
```

**After:**
```
### Quick Start: Calculate Motor Angles for a Crystal Reflection

[Show fourcv demo: define geometry, orient crystal, get angles]

### What Just Happened?

This demo used:
- **Stage definitions** to describe the diffractometer
- **Crystallographic calculations** to orient the crystal  
- **Forward diffraction** to compute motor angles

### Ready for More? Try These Demonstrations:

- **Six-circle geometry** — azimuthal control and more reciprocal-space access
- **Kappa geometry** — space-constrained setups
- **Reciprocal-space planning** — trajectory calculations for complex measurements
- [More demos...]

### Advanced Features

[Details on diffraction modes, all calculation methods, etc.]
```

---

#### **12. Establish a Demonstration Maintenance Standard**

Demonstrations are now **first-class artifacts**. Treat them accordingly:

- [ ] Every demonstration **must run without errors** (test them automatically)
- [ ] Every demonstration **must show one coherent idea** (not five features at once)
- [ ] Every demonstration **must be annotated** (explain the *why*, not just the *what*)
- [ ] Every demonstration **must be current** (update with API changes immediately)
- [ ] Every demonstration **must answer a real question** ("When would I use this?")

---

### Summary: The Reframe

**Old model:** "Here's a powerful, flexible package. Here are some example geometries to show you what's possible."

**New model:** "Here are the standard geometries you'll use. Each one demonstrates a core feature and answers a real experimental question. Together, they show you the full power of the package and how to extend it yourself."

Demonstrations become the **primary interface** to the package, not the decoration around it.