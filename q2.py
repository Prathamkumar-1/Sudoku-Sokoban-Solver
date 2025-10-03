"""
Sokoban Solver using SAT (Boilerplate)
--------------------------------------
Instructions:
- Implement encoding of Sokoban into CNF.
- Use PySAT to solve the CNF and extract moves.
- Ensure constraints for player movement, box pushes, and goal conditions.

Grid Encoding:
- 'P' = Player
- 'B' = Box
- 'G' = Goal
- '#' = Wall
- '.' = Empty space
"""

from pysat.formula import CNF
from pysat.solvers import Solver

# Directions for movement
DIRS = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

# action indices (NOOP included)
DIR_INDEX = {'U': 0, 'D': 1, 'L': 2, 'R': 3, 'NOOP': 4} #NOOP is for no operation


class SokobanEncoder:
    def _init_(self, grid, T):
        """
        Initialize encoder with grid and time limit.

        Args:
            grid (list[list[str]]): Sokoban grid.
            T (int): Max number of steps allowed.
        """
        self.grid = grid
        self.T = T
        self.N = len(grid)
        self.M = len(grid[0])

        self.goals = []
        self.boxes = []
        self.player_start = None

        # TODO: Parse grid to fill self.goals, self.boxes, self.player_start
        self._parse_grid()

        self.num_boxes = len(self.boxes)
        self.cnf = CNF()

        # variable layout sizes (base offsets)
        self.player_var_count = (T + 1) * self.num_cells
        self.box_var_count = self.num_boxes * (T + 1) * self.num_cells
        self.action_var_count = T * 5  # U,D,L,R,NOOP

        self.base_player = 1
        self.base_box = self.base_player + self.player_var_count
        self.base_action = self.base_box + self.box_var_count

        # aux variables will start after fixed vars
        self.next_aux = self.base_action + self.action_var_count
        self.max_var = self.next_aux - 1

    def _parse_grid(self):
        """Parse grid to find player, boxes, and goals."""
        # TODO: Implement parsing logic
        self.walls = set()
        for i in range(self.N):
            for j in range(self.M):
                c = self.grid[i][j]
                if c == '#':
                    self.walls.add((i, j))
                elif c == 'P':
                    self.player_start = (i, j)
                elif c == 'B':
                    self.boxes.append((i, j))
                elif c == 'G':
                    self.goals.append((i, j))
                elif c == '.':
                    pass
                else:
                    # treat unknown as floor
                    pass
        # free cells for enumerations
        self.free_cells = [(i, j) for i in range(self.N) for j in range(self.M) if (i, j) not in self.walls]
        self.num_cells = self.N * self.M

    # ---------- helper utilities ----------
    def new_aux(self):
        v = self.next_aux
        self.next_aux += 1
        if v > self.max_var:
            self.max_var = v
        return v

    # def in_bounds(self, x, y):
    #     return 0 <= x < self.N and 0 <= y < self.M

    def in_bounds_free(self, x, y):
        return 0 <= x < self.N and 0 <= y < self.M and (x, y) not in self.walls

    def at_least_one(self, lits):
        if lits:
            self.cnf.append(lits)

    def at_most_one(self, lits):
        for i in range(len(lits)):
            for j in range(i + 1, len(lits)):
                self.cnf.append([-lits[i], -lits[j]])

    def exactly_one(self, lits):
        self.at_least_one(lits)
        self.at_most_one(lits)

    def _cell_indexer(self, x, y):
        return x * self.M + y

    # ---------------- Variable Encoding ----------------
    def var_player(self, x, y, t):
        """
        Variable ID for player at (x, y) at time t.
        """
        # TODO: Implement encoding scheme
        return self.base_player + t * self.num_cells + self._cell_indexer(x, y)

    def var_box(self, b, x, y, t):
        """
        Variable ID for box b at (x, y) at time t.
        """
        # TODO: Implement encoding scheme
        stride = (self.T + 1) * self.num_cells
        return self.base_box + b * stride + t * self.num_cells + self._cell_indexer(x, y)

    def var_action(self, d, t):
        return self.base_action + t * 5 + DIR_INDEX[d]

    # ---------------- Encoding Logic ----------------
    def encode(self):
        # 1. Initial conditions
        if self.player_start is None:
            # no player -> UNSAT
            self.cnf.append([])
        else:
            self.cnf.append([self.var_player(self.player_start[0], self.player_start[1], 0)])

        for b_idx, (bx, by) in enumerate(self.boxes):
            self.cnf.append([self.var_box(b_idx, bx, by, 0)])

        # 2. Player movement
        # player exactly one cell per time
        for t in range(self.T + 1):
            self.exactly_one([self.var_player(x, y, t) for (x, y) in self.free_cells])

        # 3. Box movement (push rules)
        # each box exactly one cell per time
        for b in range(self.num_boxes):
            for t in range(self.T + 1):
                self.exactly_one([self.var_box(b, x, y, t) for (x, y) in self.free_cells])

        # 4. Non-overlap constraints
        # boxes non-overlap per cell per time
        for t in range(self.T + 1):
            for (x, y) in self.free_cells:
                self.at_most_one([self.var_box(b, x, y, t) for b in range(self.num_boxes)])

        # player and box cannot share a cell
        for t in range(self.T + 1):
            for (x, y) in self.free_cells:
                p = self.var_player(x, y, t)
                for b in range(self.num_boxes):
                    self.cnf.append([-p, -self.var_box(b, x, y, t)])

        # forbid player/boxes on walls at all times
        for t in range(self.T + 1):
            for (x, y) in self.walls:
                self.cnf.append([-self.var_player(x, y, t)])
                for b in range(self.num_boxes):
                    self.cnf.append([-self.var_box(b, x, y, t)])

        # 5) actions: exactly one action per step (U, D, L, R, NOOP)
        for t in range(self.T):
            self.exactly_one([self.var_action(d, t) for d in ['U', 'D', 'L', 'R', 'NOOP']])

        # -------------------------------------------------------
        # 7) Build transition rules for each time step
        # -------------------------------------------------------
        for t in range(self.T):   # loop over time steps
            # -----------------------------------
            # A) PLAYER MOVEMENT
            # -----------------------------------
            possible_moves_into = {cell: [] for cell in self.free_cells}  
            # (cell -> list of auxiliary vars that justify being there)

            push_records = []  # keep track of pushes for box logic later

            # For every cell where the player could be at time t+1
            for (x_next, y_next) in self.free_cells:

                # (1) NOOP: Player stays still
                aux_noop = self.new_aux()
                # If player at (x_next,y_next) at t AND action=NOOP -> aux_noop
                self.cnf.append([-self.var_player(x_next, y_next, t), -self.var_action('NOOP', t), aux_noop])
                # aux_noop -> (player was at (x_next,y_next,t), action=NOOP, and stays at (x_next,y_next,t+1))
                self.cnf.append([-aux_noop, self.var_player(x_next, y_next, t)])
                self.cnf.append([-aux_noop, self.var_action('NOOP', t)])
                self.cnf.append([-aux_noop, self.var_player(x_next, y_next, t+1)])

                possible_moves_into[(x_next, y_next)].append(aux_noop)

                # (2) WALK / PUSH from neighboring cells into (x_next,y_next)
                for direction, (dx, dy) in DIRS.items():
                    x_prev, y_prev = x_next - dx, y_next - dy  # where player came from
                    if not self.in_bounds_free(x_prev, y_prev):
                        continue

                    # (2a) WALK: player walks into (x_next,y_next) if empty
                    aux_walk = self.new_aux()
                    box_vars_here = [self.var_box(b, x_next, y_next, t) for b in range(self.num_boxes)]
                    self.cnf.append([-self.var_player(x_prev, y_prev, t), -self.var_action(direction, t)] + box_vars_here + [aux_walk])
                    # aux_walk -> player was at (x_prev,y_prev), action=direction, no box at (x_next,y_next), player ends at (x_next,y_next)
                    self.cnf.append([-aux_walk, self.var_player(x_prev, y_prev, t)])
                    self.cnf.append([-aux_walk, self.var_action(direction, t)])
                    for b in range(self.num_boxes):
                        self.cnf.append([-aux_walk, -self.var_box(b, x_next, y_next, t)])
                    self.cnf.append([-aux_walk, self.var_player(x_next, y_next, t+1)])

                    possible_moves_into[(x_next, y_next)].append(aux_walk)

                    # (2b) PUSH: player pushes a box from (x_next,y_next) to (x_next+dx,y_next+dy)
                    x_box_new, y_box_new = x_next + dx, y_next + dy
                    if not self.in_bounds_free(x_box_new, y_box_new):
                        continue

                    for b_idx in range(self.num_boxes):
                        aux_push = self.new_aux()
                        # Conditions: player at (x_prev,y_prev), action=direction, box b at (x_next,y_next),
                        # and no box blocking the next cell (x_box_new,y_box_new).
                        box_in_next_vars = [self.var_box(bb, x_box_new, y_box_new, t) for bb in range(self.num_boxes)]
                        self.cnf.append(
                            [-self.var_player(x_prev, y_prev, t),
                            -self.var_action(direction, t),
                            -self.var_box(b_idx, x_next, y_next, t)] + box_in_next_vars + [aux_push]
                        )
                        # aux_push -> all those conditions, plus effects
                        self.cnf.append([-aux_push, self.var_player(x_prev, y_prev, t)])
                        self.cnf.append([-aux_push, self.var_action(direction, t)])
                        self.cnf.append([-aux_push, self.var_box(b_idx, x_next, y_next, t)])
                        for bb in range(self.num_boxes):
                            self.cnf.append([-aux_push, -self.var_box(bb, x_box_new, y_box_new, t)])
                        self.cnf.append([-aux_push, self.var_player(x_next, y_next, t+1)])   # player steps into old box cell
                        self.cnf.append([-aux_push, self.var_box(b_idx, x_box_new, y_box_new, t+1)]) # box moves forward

                        push_records.append({'aux': aux_push, 'box': b_idx, 'from': (x_next, y_next), 'to': (x_box_new, y_box_new), 'time': t})
                        possible_moves_into[(x_next, y_next)].append(aux_push)

            # Backward link: If player at (x,y) at t+1, then some valid transition fired
            for (x, y), aux_list in possible_moves_into.items():
                if aux_list:
                    self.cnf.append([-self.var_player(x, y, t+1)] + aux_list)
                else:
                    self.cnf.append([-self.var_player(x, y, t+1)])

            # -----------------------------------
            # B) BOX MOVEMENT
            # -----------------------------------
            push_into = {}   # (box,cell) -> aux vars that push box into this cell
            push_outof = {}  # (box,cell) -> aux vars that push box away from this cell
            for info in push_records:
                push_into.setdefault((info['box'], info['to']), []).append(info['aux'])
                push_outof.setdefault((info['box'], info['from']), []).append(info['aux'])

            for b in range(self.num_boxes):
                for (qx, qy) in self.free_cells:
                    in_list  = push_into.get((b, (qx, qy)), [])
                    out_list = push_outof.get((b, (qx, qy)), [])

                    # "stay" case: box stays put if it was there and not pushed away
                    aux_stay = self.new_aux()
                    self.cnf.append([-self.var_box(b, qx, qy, t)] + out_list + [aux_stay])
                    self.cnf.append([-aux_stay, self.var_box(b, qx, qy, t)])
                    for push_out in out_list:
                        self.cnf.append([-aux_stay, -push_out])
                    self.cnf.append([-aux_stay, self.var_box(b, qx, qy, t+1)])

                    # Backward: if box is at (qx,qy) at t+1, then either it stayed or was pushed in
                    ors = [aux_stay] + in_list
                    if ors:
                        self.cnf.append([-self.var_box(b, qx, qy, t+1)] + ors)
                    else:
                        self.cnf.append([-self.var_box(b, qx, qy, t+1)])


        # 8) goal condition at time T: every box must be on some goal
        if self.num_boxes > 0 and len(self.goals) == 0:
            # impossible
            self.cnf.append([])

        for b_idx in range(self.num_boxes):
            glits = [self.var_box(b_idx, gx, gy, self.T) for (gx, gy) in self.goals]
            if glits:
                self.at_least_one(glits)
            else:
                self.cnf.append([])

        return self.cnf


def decode(model, encoder):
    """
    Decode SAT model into list of moves ('U', 'D', 'L', 'R').

    Args:
        model (list[int]): Satisfying assignment from SAT solver.
        encoder (SokobanEncoder): Encoder object with grid info.

    Returns:
        list[str]: Sequence of moves.
    """
    N, M, T = encoder.N, encoder.M, encoder.T

    # TODO: Map player positions at each timestep to movement directions
    true_set = set(l for l in model if l > 0)
    moves = []
    for t in range(T):
        found = False
        # exactly one action was enforced, so choosing it
        for d in ['U', 'D', 'L', 'R', 'NOOP']:
            if encoder.var_action(d, t) in true_set:
                if d != 'NOOP':
                    moves.append(d)
                found = True
                break
        if not found:
            # broken model / decode fail
            return -1
    return moves  # return list as tester expects


# ---------- public solve function ----------
def solve_sokoban(grid, T):
    """
    DO NOT MODIFY THIS FUNCTION.

    Solve Sokoban using SAT encoding.

    Args:
        grid (list[list[str]]): Sokoban grid.
        T (int): Max number of steps allowed.

    Returns:
        list[str] or "unsat": Move sequence or unsatisfiable.
    """
    encoder = SokobanEncoder(grid, T)
    cnf = encoder.encode()

    with Solver(name='g3') as solver:
        solver.append_formula(cnf)
        if not solver.solve():
            return -1

        model = solver.get_model()
        if not model:
            return -1

        return decode(model, encoder)