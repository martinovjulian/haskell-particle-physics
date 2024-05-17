module Main (main) where

import Linear.Metric
import Linear.V3
import Linear.Vector
import Control.Parallel.Strategies
import Vis


{- grid size, i.e. initializes n x n x n cube of atoms (don't use n>5 for animated executions; for n=5, you will likely need to use N=4 cores)-}
n :: Int
n = 4

 {- chunk size to use in parListChunk. Set this to numCores / n^3 -}
chunkSize :: Int
chunkSize = 32

{- cubical container side length; if you have a large n, e.g. >=9, make sure to expand the box so there is space for the particles to take discrete time steps without getting unphysically close together, e.g. s >= 10 -}
s :: Float
s = 10
 

main :: IO ()
main = mainAnim {-choose from mainAnim or mainNoAnim to either run the
                    simulation 3D animated in a GUI window, or run a finite number
                    of time steps without any animation-}

{- more parameters; you probably don't want to change these -}
rad :: Float
rad = 0.15 {- atom radius -}
timeStep :: Float
timeStep = 0.1 {- time step length -}

{- define the Atom data type, the basic unit of our simulation -}
data Atom = Atom { i :: Int,       -- index in the array
                   r :: V3 Float,  -- position vector
                   _v :: V3 Float } -- velocity vector


{- position update in linear time -}
rstep :: Float -> Atom -> V3 Float -> Atom
rstep timeStepR (Atom i r v) a = Atom i r' v
  where r' = r ^+^ (v ^* timeStepR) ^+^ (0.5*timeStepR**2 *^ a)

{- velocity update in linear time -}
vstep :: Float -> Atom -> V3 Float -> Atom
vstep timeStepV atom a = Atom i r v'
  where (Atom i r v) = atom
        v' = (bound atom) * (v + (0.5 * timeStepV) *^ a)
        bound (Atom _ (V3 x y z) _) = V3 xf yf zf -- enforces rigid wall boundary condition
          where xf = if (abs x + rad > s/2) then (-1) else 1
                yf = if (abs y + rad > s/2) then (-1) else 1
                zf = if (abs z + rad > s/2) then (-1) else 1

{- force update in quadratic time -}
fstep :: [Atom] -> [V3 Float]
fstep atoms = fTot atoms atoms
  where fTot [a] bs = fOne a bs
        fTot (a:as) bs = fOne a bs ^+^ fTot as bs
        fTot [] _ = fTot atoms atoms

{- helper function to calculate the total net force acting on a single atom -}
fOne :: Atom -> [Atom] -> [V3 Float]
fOne atom = fmap (f atom)
  where f a b = if (i a==i b) then 0 else lennardJones a b
        lennardJones a b = (1 / (norm d)^(14 :: Integer) - 0.5 / (norm d)^(8 :: Integer)) *^ d
          where d = (r b) ^-^ (r a)

{- velocity Verlet algorithm -}
step :: Float -> [Atom] -> [Atom]
step timeStepS atoms = zipWith (vstep timeStepS) r' f'
  where f = fstep atoms `using` parListChunk chunkSize rdeepseq
        r' = zipWith (rstep timeStepS) atoms f
        f' = f ^+^ (fstep r' `using` parListChunk chunkSize rdeepseq)

{- run the program with animation enabled -}
mainAnim :: IO ()
mainAnim = simulate options refreshRate initConfig draw update
  where options =
          ( defaultOpts
            { optWindowName = "Lennard Particle Sim",
              optBackgroundColor = Just white,
              optWindowSize = Just (1280, 720),
              optWindowPosition = Just (160, 50),
              optAntialiasing = Smoothed
            }
          )
        refreshRate = 0.01
        initConfig = grid n
        draw config = VisObjects $ [box] ++ (drawAtom <$> config `using` parListChunk chunkSize rseq)
          where box = Trans (V3 0 0 0) $ Box (s, s, s) Wireframe black
                drawAtom atom = Trans (r atom) $ Sphere rad Solid aquamarine
        update _ config = step timeStep config


{- initialize a n x n x n cubical grid as the initial atom configuration -}
grid :: Int -> [Atom]
grid nGrid = zipWith3 Atom [1..(nGrid^(3 :: Integer))] (cube nGrid nGrid) (replicate (nGrid^(3 :: Integer)) (V3 0 0 0))
  where cube _ 0 = []
        cube d i = square d d z ++ cube d (i-1)
          where z = s/2 - (fromIntegral i * s/fromIntegral (d+1))

square :: Int -> Int -> Float -> [V3 Float]
square _ 0 _ = []
square d i z = row d d y ++ square d (i-1) z
  where y = s/2 - (fromIntegral i * s/fromIntegral (d+1))
        row _ 0 _ = []
        row dSquare iSquare ySquare = V3 x y z : row dSquare (iSquare-1) ySquare
          where x = s/2 - (fromIntegral iSquare * s/fromIntegral (dSquare+1))