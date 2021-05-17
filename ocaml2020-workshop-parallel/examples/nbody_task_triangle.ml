
module T = Domainslib.Task

type planet = { mutable x : float;  mutable y : float;  mutable z : float;
                mutable vx: float;  mutable vy: float;  mutable vz: float;
                mass : float }

(* Compact index into flat array (of size) for
   matrix of upper triangle (including diagonal)
   0th row contains n elements,
   1st row contains n-1 elements,
   kth row contains n-k elements
*)
let index n i j = (i * ( (2*n) - 1 - i) / 2) + j

let initialize_velocities n_bodies =
  Array.init (index n_bodies n_bodies n_bodies) (fun _ -> {x=0.; y=0.; z=0.; vx=0.; vy=0.; vz=0.; mass=0.})

let advance pool bodies velocities dt =
  let n_bodies = Array.length bodies in
  (* calculate velocity increments *)
  T.parallel_for pool
    ~chunk_size:1
    ~start:0
    ~finish:(n_bodies - 1)
    ~body:(fun i ->
      let b = bodies.(i) in
      let vx = ref 0. and vy = ref 0. and vz = ref 0. in
      for j = i+1 to n_bodies - 1 do
        let b' = bodies.(j) in
        let dx = b.x -. b'.x and dy = b.y -. b'.y and dz = b.z -. b'.z in
        let dist2 = dx *. dx +. dy *. dy +. dz *. dz in
        let mag = dt /. (dist2 *. sqrt(dist2)) in

        vx := !vx -. dx *. b'.mass *. mag;
        vy := !vy -. dy *. b'.mass *. mag;
        vz := !vz -. dz *. b'.mass *. mag;

        let vb' = velocities.(index n_bodies i j) in
        vb'.vx <- dx *. b.mass *. mag;
        vb'.vy <- dy *. b.mass *. mag;
        vb'.vz <- dz *. b.mass *. mag
      done;
      let vb = velocities.(index n_bodies i i) in
      vb.vx <- !vx; vb.vy <- !vy; vb.vz <- !vz
    );
  let k = ref 0 in
  (* accumulate velocities *)
  for i = 0 to n_bodies - 1 do
    for j = i to n_bodies - 1 do
      let b = bodies.(j) and v = velocities.(!k) in
      b.vx <- b.vx +. v.vx;
      b.vy <- b.vy +. v.vy;
      b.vz <- b.vz +. v.vz;
      incr k
    done
  done;
  (* advance positions *)
  for i = 0 to n_bodies - 1 do
    let b = bodies.(i) in
    b.x <- b.x +. dt *. b.vx;
    b.y <- b.y +. dt *. b.vy;
    b.z <- b.z +. dt *. b.vz;
  done

let energy bodies =
  let e = ref 0. in
  for i = 0 to Array.length bodies - 1 do
    let b = bodies.(i) in
    e := !e +. 0.5 *. b.mass *. (b.vx *. b.vx +. b.vy *. b.vy +. b.vz *. b.vz);
    for j = i+1 to Array.length bodies - 1 do
      let b' = bodies.(j) in
      let dx = b.x -. b'.x  and dy = b.y -. b'.y  and dz = b.z -. b'.z in
      let distance = sqrt(dx *. dx +. dy *. dy +. dz *. dz) in
      e := !e -. (b.mass *. b'.mass) /. distance
    done
  done;
  !e

let pi = 3.141592653589793
let solar_mass = 4. *. pi *. pi
let days_per_year = 365.24

let offset_momentum bodies =
  let px = ref 0. and py = ref 0. and pz = ref 0. in
  for i = 0 to Array.length bodies - 1 do
    px := !px +. bodies.(i).vx *. bodies.(i).mass;
    py := !py +. bodies.(i).vy *. bodies.(i).mass;
    pz := !pz +. bodies.(i).vz *. bodies.(i).mass;
  done;
  bodies.(0).vx <- -. !px /. solar_mass;
  bodies.(0).vy <- -. !py /. solar_mass;
  bodies.(0).vz <- -. !pz /. solar_mass

let initialize_bodies n_bodies =
  Array.init n_bodies (fun _ ->
    { x = (Random.float  10.);
      y = (Random.float 10.);
      z = (Random.float 10.);
      vx= (Random.float 5.) *. days_per_year;
      vy= (Random.float 4.) *. days_per_year;
      vz= (Random.float 5.) *. days_per_year;
      mass=(Random.float 10.) *. solar_mass; })

let () =
  let n = int_of_string(Sys.argv.(1)) in
  let n_bodies = int_of_string(Sys.argv.(2)) in
  let n_domains = int_of_string(Sys.argv.(3)) in
  let pool = T.setup_pool ~num_domains:(n_domains - 1) in
  let bodies = initialize_bodies n_bodies in
  offset_momentum bodies;
  Printf.printf "%.9f\n" (energy bodies);
  let velocities = initialize_velocities n_bodies in
  for _i = 1 to n do advance pool bodies velocities 0.01 done;
  Printf.printf "%.9f\n" (energy bodies);
  T.teardown_pool pool
