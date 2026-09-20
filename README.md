# quakem
robust movement systems, and i demonstrate two of them which approach the problem differently. using the two systems just means switching the collider type to either cylinder or capsule.
these systems should be very easy to embed beyond just quake movement, such as server authority or a base for custom player physics. these approaches are much more deterministic.

## cylinder
* uses a probably novel separating axis theorem extension that adds a parabolic TOI (works best on exact primitives, can run faster than FEV/VClip in this scenario)
* uses a greedy active set solver for discrete collision detection
* analytical cylinder tracing, supports tracing with velocity and acceleration
* physics interactions supported, although main purpose is for accuracy
* does not include terrain/meshparts because i couldn't really find any clean exact solutions that mesh well with roblox's api
* best for more competitive movement (bhop, surf)
* hopefully no rampbugs
* framerate independent

## capsule
* supports terrain + meshparts (concave approximation through planes)
* approximates collision through a linear sweep, but uses the greedy active set solver
* physics interactions are fully supported
* best for general environments
* hopefully no rampbugs
* framerate independent

## demo
[surf utopia](https://www.youtube.com/watch?v=5fJUCuLyTMk&feature=youtu.be)