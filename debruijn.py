import sympy
import itertools
import time
from sympy.abc import x


def powerset(iterable, reverse=True):
    """powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"""
    s = list(iterable)
    if not reverse:
        return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))
    else:
        return itertools.chain.from_iterable(itertools.combinations(s, r) for r in reversed(range(len(s)+1)))


def hamming_weight(n):
    c = 0
    while n:
        c += 1
        n &= n - 1
    return c


def is_primitive(poly):
    if not poly.is_irreducible:
        return False

    degree = int(sympy.degree(poly))
    if sympy.isprime(2**degree - 1):
        return True

    # compare with generating the full sequence
    for k in (d for d in sympy.divisors(2**degree-1) if degree < d < (2**degree - 1)):
        q = sympy.Poly(x**k + 1, x, modulus=2)
        if sympy.rem(q, poly, modulus=2).is_zero:
            return False
    return True


def generate_primitives(degree):
    if degree == 1:
        return [sympy.Poly(x+1, x, modulus=2)]

    for k in xrange(1, 2**(degree-1)):
        if hamming_weight(k) % 2 == 1:
            if int(bin(k)[:1:-1] + '0' * (degree - int(sympy.log(k, 2)) - 2), 2) < k:
                continue
            proto_poly = 0
            power = 0
            while k:
                if k % 2:
                    proto_poly += x**power
                k >>= 1
                power += 1
            proto_poly *= x
            proto_poly += x**degree + 1
            poly = sympy.Poly(proto_poly, x, modulus=2)
            if is_primitive(poly):
                break

    decimations = []
    for k in range(1, 2**(degree-1)):
        if sympy.gcd(k, 2**degree - 1) != 1:
            continue
        is_new_coset = True
        for i in range(1, degree):
            if (k * 2**i) % (2**degree - 1) in decimations:
                is_new_coset = False
                break
        if is_new_coset:
            decimations.append(k)

    return sorted(map(lambda a: poly_decimation(poly, a), decimations), key=lambda a: a.all_coeffs())


def get_associate_poly(poly_list):
    associates = []
    indexes = []
    for d in set(map(sympy.degree, poly_list)):
        poly_subset = [poly for poly in poly_list if sympy.degree(poly) == d]
        primitive_polys = generate_primitives(d)
        for poly in poly_subset:
            indexes.append(poly_list.index(poly))
            if poly in primitive_polys:
                associates.append({'poly': poly,
                                   'order': 1,
                                   'period': 2**d-1,
                                   'associate': poly})
                continue

            # routine to find sequence's period
            sequence = [1] + [0]*(d-1)
            e = 1
            divisor = sympy.divisors(2**d-1)
            for e in divisor:
                if e < d:
                    continue
                while len(sequence) < 2*e:
                    sequence.append(lfsr_from_poly(poly, sequence[-d:])[-1])

                is_periodic = True
                for i in range(e):
                    is_periodic = is_periodic and sequence[i] == sequence[i+e]
                    if not is_periodic:
                        is_periodic = False
                        break
                if is_periodic:
                    break

            # routine to find associate polynomial
            for primitive in primitive_polys:
                if poly_decimation(primitive, (2**d-1)/e) == poly:
                    associates.append({'poly': poly,
                                       'order': (2**d-1)/e,
                                       'period': e,
                                       'associate': primitive})
                    break
    return [p for _, p in sorted(zip(indexes, associates))]


# THIS SECTION IS LIFTED FROM ZECH.PY #
def lfsr_from_poly(poly, state):
    lfsr = list(reversed(poly.all_coeffs()))[:-1]
    next_state = state + [sum(map(lambda a, b: a & int(b), state, lfsr)) % 2]
    return next_state[1:]


def seq_decimation(p, t, offset=0, c_state=None):
    deg = sympy.degree(p)
    if c_state is None:
        c_state = [1] * deg
    ret = [0] * (2*deg)
    for _ in range(offset):
        c_state = lfsr_from_poly(p, c_state)
    indexes = map(lambda a: (a*t) % (2**deg - 1), range(2*deg))
    ctr = 0
    i = 0
    while ctr < 2*deg:
        if i in indexes:
            pos = indexes.index(i)
            ret[pos] = c_state[0]
            indexes[pos] = -1
            ctr += 1
        c_state = lfsr_from_poly(p, c_state)
        i += 1
        if i >= 2**deg - 1:
            i -= 2**deg - 1
    return ret


def poly_decimation(p, t):
    from operator import mul
    n = sympy.degree(p)
    s = seq_decimation(p, t)
    while len(s) < 2*n:
        s += s

    cd = sympy.Poly(1, x, modulus=2)
    l, m, bd = 0, -1, 1
    for i in range(2*n):
        sub_cd = list(reversed(cd.all_coeffs()))[1:l+1]
        sub_s = list(reversed(s[i-l:i]))
        sub_cd += [0] * (len(sub_s) - len(sub_cd))
        disc = s[i] + sum(map(mul, sub_cd, sub_s))
        if disc % 2 == 1:
            td = cd
            cd += bd * sympy.Poly(x**(i - m), x, modulus=2)
            if l <= i/2:
                l = i + 1 - l
                m = i
                bd = td
    if sympy.degree(cd) == n:
        cd = sympy.Poly(reversed(cd.all_coeffs()), x, modulus=2)
        return cd
    else:
        return None


def get_spanning_trees(mat):
    """Initialize the spanning tree generator and return it."""
    d1_index = [i for i in range(mat.rows) if mat[i, i] == 1]
    if d1_index:
        init_edges = []
        init_verts = []
        for index in d1_index:
            if index not in init_verts:
                n_index = list(mat[:, index]).index(-1)
                init_edges.append((index, n_index))
                init_verts += [n_index]
        else:
            init_verts.sort()
        init_verts = d1_index + init_verts[:1]
        return span_tree_gen(mat, init_verts, init_edges, len(d1_index))
    else:
        init_verts = [0]
        return span_tree_gen(mat, init_verts)


def span_tree_gen(mat, verts_in, edgelist=None, cur_vert=0):
    """Generate and yield all spanning trees of a given graph."""
    if len(edgelist) == mat.rows-1:
        yield sorted(edgelist)
    else:
        if edgelist is None:
            edgelist = []
        if cur_vert < len(verts_in):
            neighbor = mat[verts_in[cur_vert], :]
            nlist = [i for i, j in enumerate(neighbor) if j == -1 if i not in verts_in]
            for vertices in powerset(nlist):
                verts = list(vertices)
                new_edges = [(verts_in[cur_vert], i) for i in verts]
                for next_tree in span_tree_gen(mat, verts_in + verts, edgelist + new_edges, cur_vert + 1):
                    yield next_tree
# END SECTION LIFTED FROM ZECH.PY #


class DeBruijnPoly:
    def __init__(self, *args):
        """ takes binary strings which are converted into polys """
        if not args:
            raise ValueError('no arguments passed (at least 1 expected)')

        # properties modifiable by user
        self.single_mode = False
        self.log = False
        self._initial_state = None

        # properties that are read-only
        self._initialized = False
        self._bit_exhausted = False
        self._polys = []
        self._states = []
        self._associates = []
        self._adjacency_dict = {}
        self._ordered_params = []
        self._degree = 0
        self._poly = None
        self._p_matrix = None
        self._adjacency_matrix = None
        self._param_generator = None
        self._db_seq_param = None
        self._anf = None
        self._anf_degree = 0

        for binary_string in args:
            binary_seq = map(lambda a: 0 if a == '0' else 1, binary_string)
            proto_poly = reduce(lambda a, b: a*x + b, binary_seq, 0)
            poly = sympy.Poly(proto_poly, x, modulus=2)
            if poly.is_irreducible and poly not in self._polys:
                self._polys.append(poly)
        self._polys.sort(key=sympy.degree)
        self._poly = reduce(lambda a, b: a * b, self._polys)

        self._degree = sum(map(sympy.degree, self._polys))
        self._initial_state = [0] * self._degree

        self._sym = [sympy.Symbol('x_{}'.format(i), integer=True) for i in range(self._degree)]

    def initialize(self):
        """ takes a seeded DeBruijn object and prepares it for generating sequences """
        if self._initialized:
            return

        # populate states
        if self.log:
            start = time.clock()
        self._associates = get_associate_poly(self._polys)
        for entry in self._associates:
            entry_state = []
            degree = sympy.degree(entry['associate'])
            init_state = [1] * degree
            for i in xrange(entry['order']):
                entry_state.append(seq_decimation(entry['associate'], entry['order'], i, init_state)[:degree])
            entry_state.append([0] * degree)
            self._states.append(entry_state)
        if self.log:
            print 'associates and states : {} s'.format(time.clock() - start)

        # find special state
        if self.log:
            start = time.clock()
        p_matrix = []
        for poly in self._polys:
            degree = sympy.degree(poly)
            for i in xrange(degree):
                state = [0] * degree
                state[i] = 1
                for j in xrange(self._degree - degree):
                    state.append(lfsr_from_poly(poly, state[-degree:])[-1])
                p_matrix += state
        p_matrix = sympy.Matrix(self._degree, self._degree, p_matrix)
        self._p_matrix = p_matrix
        special_state = map(int, sympy.Matrix(1, self._degree, [1] + [0] * (self._degree - 1)) * p_matrix.inv_mod(2))
        special_states = []
        i = 0
        for poly in self._polys:
            special_states.append(special_state[i:i+sympy.degree(poly)])
            i += sympy.degree(poly)
        if self.log:
            print 'special states        : {} s'.format(time.clock() - start)

        # find viable pairs
        # notes: all_pairs = list of dictionary for each polynomial, where the keys are pairs of states,
        #                    and the entries are the corresponding shifts for the states.
        if self.log:
            start = time.clock()
        all_pairs = []
        for i, entry in enumerate(self._associates):
            cur_pairs = {}
            special_list = [special_states[i]]
            for _ in xrange(entry['period'] - 1):
                special_list.append(lfsr_from_poly(entry['poly'], special_list[-1]))
            for state_1 in xrange(entry['order']+1):
                added_state = map(lambda a: map(lambda b, c: b ^ c, self._states[i][state_1], a), special_list)
                for state_2 in xrange(state_1, entry['order']+1):
                    cur_state = self._states[i][state_2][:]
                    for shift_2 in xrange(entry['period'] if state_2 != entry['order'] else 1):
                        if cur_state in added_state:
                            shift_1 = added_state.index(cur_state)
                            if state_2 == entry['order']:
                                shift_2 = shift_1
                            if (state_1, state_2) in cur_pairs:
                                cur_pairs[(state_1, state_2)].append((-shift_1 % entry['period'],
                                                                      (shift_2 - shift_1) % entry['period']))
                                # either include flips here, or in conjugate generator
                                if state_1 != state_2:
                                    cur_pairs[(state_2, state_1)].append(((shift_2 - shift_1) % entry['period'],
                                                                          -shift_1 % entry['period']))
                            else:
                                cur_pairs[(state_1, state_2)] = [(-shift_1 % entry['period'],
                                                                  (shift_2 - shift_1) % entry['period'])]
                                if state_1 != state_2:
                                    cur_pairs[(state_2, state_1)] = [(((shift_2 - shift_1) % entry['period'],
                                                                       -shift_1 % entry['period']))]
                            cur_state = lfsr_from_poly(entry['poly'], cur_state)
                        else:
                            cur_state = lfsr_from_poly(entry['poly'], cur_state)
            all_pairs.append(cur_pairs)
        if self.log:
            print 'find pairs            : {} s'.format(time.clock() - start)

        # find conjugate pairs and construct adjacency graph
        if self.single_mode:
            param_list = []
            for n in xrange(2**len(self._polys)):
                lookup = [n & 2**i for i in range(len(self._polys))]
                periods = [self._associates[i]['period'] if lookup[i] else 1 for i in range(len(self._polys))]
                orders = [self._associates[i]['order'] if lookup[i] else 1 for i in range(len(self._polys))]
                products = []
                for i in xrange(len(self._polys)):
                    products.append(sympy.gcd(periods[i], reduce(sympy.lcm, periods[:i], 1)))
                ranges = [range(orders[i]) if lookup[i] else [self._associates[i]['order']]
                          for i in range(len(self._polys))]
                param_list += list(itertools.product(*(ranges + map(range, products))))
            param_list.sort(reverse=True)
            self._ordered_params = param_list
            connection = [[i] for i in range(len(param_list))]

        if self.log:
            start = time.clock()
        adjacency_dict = self._adjacency_dict
        for param_1, param_2, shifts_1, shifts_2 in self.__conjugate_pair_generator(all_pairs):
            if param_1 != param_2:
                if param_1 not in adjacency_dict:
                    adjacency_dict[param_1] = {}
                if param_2 not in adjacency_dict[param_1]:
                    adjacency_dict[param_1][param_2] = []

                if not self.single_mode:
                    adjacency_dict[param_1][param_2].append(shifts_1)
                else:
                    if param_2 not in adjacency_dict:
                        adjacency_dict[param_2] = {}
                    if param_1 not in adjacency_dict[param_2]:
                        adjacency_dict[param_2][param_1] = []

                    if shifts_1 in adjacency_dict[param_1][param_2]:
                        continue

                    adjacency_dict[param_1][param_2].append(shifts_1)
                    adjacency_dict[param_2][param_1].append(shifts_2)
                    pos = [i for i, subset in enumerate(connection) if
                           param_list.index(param_1) in subset or param_list.index(param_2) in subset]
                    if len(pos) == 1:
                        continue
                    else:
                        pos_1, pos_2 = pos
                        connection[pos_1] += connection.pop(pos_2)

                    if len(connection) == 1:
                        break

        if self.log:
            print 'conjugate pairs       : {} s'.format(time.clock() - start)

        # clean up
        del all_pairs

        if self.log:
            start = time.clock()
        if not self.single_mode:
            param_list = sorted(adjacency_dict.keys(), reverse=True)
            self._ordered_params = param_list
        adjacency_matrix = sympy.zeros(len(param_list))
        for param_1 in adjacency_dict:
            for param_2 in adjacency_dict[param_1]:
                pos_1 = param_list.index(param_1)
                pos_2 = param_list.index(param_2)
                adjacency_matrix[pos_1, pos_2] = -len(adjacency_dict[param_1][param_2])
        for i in range(adjacency_matrix.rows):
            adjacency_matrix[i, i] = -sum(adjacency_matrix[:, i])
        self._adjacency_matrix = adjacency_matrix
        if self.log:
            print 'adjacency graph       : {} s'.format(time.clock() - start)

        simplified_matrix = adjacency_matrix.applyfunc(lambda a: -1 if a < 0 else a)
        for i in range(simplified_matrix.rows):
            simplified_matrix[i, i] -= sum(simplified_matrix[i, :])

        # find spanning trees
        self._param_generator = self.__param_generator(get_spanning_trees(simplified_matrix))

        self._initialized = True
        self.refresh_bits()
        if self.log:
            print

    def __conjugate_pair_generator(self, all_pairs):
        # notes: pairs_product = cartesian product of all the pairs of states for each polynomials,
        #                        i.e. it's a tuple of size s, where s is the number of polynomials
        #        pairs = a tuple of size s, entry at index i is a state-pair for i-th polynomial.
        pairs_product = itertools.product(*map(lambda a: a.keys(), all_pairs))
        for pairs in pairs_product:
            for shifts in itertools.product(*[all_pairs[i][pair] for i, pair in enumerate(pairs)]):
                # there's probably a more efficient way of doing this
                periods_1 = [entry['period'] if pairs[i][0] != entry['order'] else 1
                             for i, entry in enumerate(self._associates)]
                periods_2 = [entry['period'] if pairs[i][1] != entry['order'] else 1
                             for i, entry in enumerate(self._associates)]
                found = True
                for l1 in itertools.product(*[range(a) for a in periods_1]):
                    for l2 in itertools.product(*[range(a) for a in periods_2]):
                        found = True
                        for (i, j) in itertools.combinations(range(len(self._polys)), 2):
                            if sympy.gcd(periods_1[i], periods_2[j]) != 1:
                                diff_1 = shifts[i][0] - shifts[j][0]
                                diff_2 = shifts[i][1] - shifts[j][1]
                                found &= (diff_1 - l1[i] + l1[j]) % sympy.gcd(periods_1[i], periods_2[j]) == 0
                                found &= (diff_2 - l2[i] + l2[j]) % sympy.gcd(periods_1[i], periods_2[j]) == 0
                                if not found:
                                    break
                        if found:
                            # we can reduce number of pairs by messing around in here, possibly
                            x_list_1 = []
                            x_list_2 = []
                            for i in range(len(self._polys)):
                                x_list_1.append(sympy.gcd(periods_1[i], reduce(sympy.lcm, periods_1[:i], 1)))
                                x_list_2.append(sympy.gcd(periods_2[i], reduce(sympy.lcm, periods_2[:i], 1)))
                            eff_shift_1 = [(shifts[i][0] - shifts[0][0]) % x_list_1[i]
                                           for i in range(len(self._polys))]
                            eff_shift_2 = [(shifts[i][1] - shifts[0][1]) % x_list_2[i]
                                           for i in range(len(self._polys))]
                            param_1 = tuple([a[0] for a in pairs] + eff_shift_1)
                            param_2 = tuple([a[1] for a in pairs] + eff_shift_2)
                            yield param_1, param_2,\
                                  [shifts[i][0] for i in range(len(self._polys))],\
                                  [shifts[i][1] for i in range(len(self._polys))]
                            break
                    if found:
                        break

    def __param_generator(self, trees):
        for tree in trees:
            edge_weights = [len(self._adjacency_dict[self._ordered_params[i]][self._ordered_params[j]])
                            for (i, j) in tree]
            weight_list = [range(w) for w in edge_weights]

            for weight in itertools.product(*weight_list):
                yield tree, weight

    def __get_algebraic_normal_form(self):
        terms = [a[0][0] for a in self._poly.terms() if a[0][0] != self._degree]
        anf = sum([self._sym[a] for a in terms])

        # pass the state_set into bits, probably by storing as bound variable?
        tree, weight = self._db_seq_param
        for w, (i, j) in enumerate(tree):
            param_1 = self._ordered_params[i]
            param_2 = self._ordered_params[j]
            state = []
            cur_state = list(param_1)[:len(self._polys)]
            for k in range(len(self._polys)):
                sub_state = self._states[k][cur_state[k]]
                for l in range(self._adjacency_dict[param_1][param_2][weight[w]][k]):
                    sub_state = lfsr_from_poly(self._polys[k], sub_state)
                state += sub_state
            state = (sympy.Matrix(1, self._degree, state) * self._p_matrix).applyfunc(lambda a: a % 2)[:]

            anf += reduce(lambda a, b: a * b, [self._sym[a] + state[a] + 1 for a in range(1, self._degree)])

        for subset in powerset(range(1, self._degree), reverse=False):
            anf += reduce(lambda a, b: a * b, map(lambda a: self._sym[a], list(subset)), 1)

        anf = sympy.poly(anf, modulus=2)
        self._anf = anf.as_expr()
        self._anf_degree = anf.total_degree()

    def bits(self):
        """ exposes a generator interface that emits bits of debruijn sequence """
        if not self._initialized:
            raise RuntimeError('object not initialized yet')

        if not self._bit_exhausted:
            tree, weight = self._db_seq_param
            state_set = []
            for w, (i, j) in enumerate(tree):
                param_1 = self._ordered_params[i]
                param_2 = self._ordered_params[j]
                state = []
                cur_state = list(param_1)[:len(self._polys)]
                for k in range(len(self._polys)):
                    sub_state = self._states[k][cur_state[k]]
                    for l in range(self._adjacency_dict[param_1][param_2][weight[w]][k]):
                        sub_state = lfsr_from_poly(self._polys[k], sub_state)
                    state += sub_state
                state = (sympy.Matrix(1, self._degree, state) * self._p_matrix).applyfunc(lambda a: a % 2)[:]
                state_set.append(state[1:])

            state = self._initial_state
            for _ in state:
                yield _

            for _ in xrange(2**self._degree - self._degree):
                if state[1:] in state_set:
                    state = lfsr_from_poly(self._poly, state)
                    state[-1] = 1 - state[-1]
                else:
                    state = lfsr_from_poly(self._poly, state)
                yield state[-1]

            self._bit_exhausted = True

    def refresh_bits(self):
        """ refreshes the bit generator when it is exhausted, moving on to the next sequence """
        try:
            params = self._param_generator.next()
        except StopIteration:
            self._anf = None
            self._anf_degree = float('-Inf')
            return -1
        else:
            self._db_seq_param = params
            self._bit_exhausted = False
            self.__get_algebraic_normal_form()
            return 0

    def sequences(self):
        """ exposes a generator interface that emits entire sequences """
        sequence = []
        while True:
            for bit in self.bits():
                sequence.append(bit)
            yield sequence
            sequence = []
            if self.refresh_bits() == -1:
                return

    @property
    def initialized(self):
        return self._initialized

    @property
    def ready(self):
        return self._initialized and not self._bit_exhausted

    @property
    def polys(self):
        return self._polys

    @property
    def degree(self):
        return self._degree

    @property
    def adjacency_matrix(self):
        return self._adjacency_matrix

    @property
    def initial_state(self):
        return self._initial_state

    @property
    def anf(self):
        return self._anf

    @property
    def anf_degree(self):
        return self._anf_degree

    @initial_state.setter
    def initial_state(self, iterable):
        state = map(lambda x: 1 if x else 0, iterable)
        if len(state) < self._degree:
            state += [0] * (self._degree - len(state))
        elif len(state) > self._degree:
            state = state[:self._degree]
        self._initial_state = state


if __name__ == '__main__':
    num_seq = 10
    print 'testing for generation of 10 debruijn sequences...'
    print

    ss = time.clock()
    tt = DeBruijnPoly('11111')  # ''11', '111', '11111', '1001001')
    tt.log = True
    # tt.single_mode = True
    tt.initialize()
    print 'number of cycles         : {} '.format(tt.adjacency_matrix.rows)
    dbseq = tt.adjacency_matrix[1:, 1:].det()
    print 'number of db seqs        : {}\n                         = 2^{}'.format(dbseq, sympy.log(dbseq, 2).evalf())
    sm = tt.adjacency_matrix.applyfunc(lambda a: -1 if a < 0 else a)
    for ii in range(sm.rows):
        sm[ii, ii] -= sum(sm[ii, :])
    dbseq = sm[1:, 1:].det()
    print 'number of spanning trees : {}\n                         = 2^{}'.format(dbseq, sympy.log(dbseq, 2).evalf())
    print ''
    for seq in tt.sequences():
        num_seq -= 1
        print ''.join(map(str, seq)[:80]) + ('...' if len(seq) > 80 else '')
        if num_seq == 0:
            break
    print '{} seconds elapsed.'.format(time.clock() - ss)
    print
