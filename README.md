# quakem
custom quake movement recreation using two approaches

(also sorry for my movement, i'm kinda rusty lol)

https://github.com/user-attachments/assets/2929f814-b3b7-4546-9515-4a1ce3b4a694

https://github.com/user-attachments/assets/7223df9a-1dfb-45ae-ab73-edbf20978720

## rationale
demonstrates how custom player colliders can be implemented, where i base my custom colliders off of two different approaches. they probably aren't like production ready colliders but i believe they are highly educational.
it's also to experiment player movement collision other than what's already out there. so, as of rn, this repo omits v-clip, (which I *think* roblox bhop employs a custom variant of in its face crawler for tracing parabolic paths,) the canonical quickhull + scene trace approach seen in quake derivatives + chickynoid, or the iterative computational geo algs (like the MPR + GJK trace combo for Bolt's move and slide.)

## features
* both colliders use custom physics, so no existing physics api (other than querying for primitives) really needs to exist for either of these to work
* repo contains two robust movement systems, and i demonstrate two different approaches.
    * capsule is a more general/approximate solution
    * cylinder is focused on primitives and analytical solutions for TOI and contacts.
* the collider options should be very easy to embed beyond just quake movement, since it can be used for server authority or as a base for custom player physics. moreover, both of these approaches are deterministic.
* designed so you can just pick and choose one of the colliders. there's some code duplication as a result (particularly with physics response) but it's just for cleaner packaging.

## cylinder
* uses an extension of conservative advancement that adds a parabolic TOI (works best on exact primitives, can run faster than V-Clip on primitives)
* uses a greedy active set solver for discrete collision detection
* analytical cylinder tracing, supports tracing with velocity and acceleration
* physics interactions are mostly supported
* strictly primitive contacts, so it does not include terrain/meshparts because i couldn't really find any clean exact solutions that mesh well with roblox's api
* best for more competitive movement (bhop, surf, kz, ...)
* hopefully no rampbugs
* framerate independent

## capsule
* supports terrain + meshparts (concave approximation through planes)
* approximates collision through a linear sweep, but uses the greedy active set solver
* physics interactions are fully supported
* best for general environments
* hopefully no rampbugs
* framerate independent

## usage
* example script in init.client.luau (uses variable dt per demonstration, you'd use fixed update hz for server auth.)
* also if you do use this, i'd appreciate it if you credit this project!! :D
