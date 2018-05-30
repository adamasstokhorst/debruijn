import debruijn as db


# create de bruijn generator object, seed it with polynomials:
#   x + 1, x^3 + x + 1, x^4 + x^3 + x^2 + x + 1.
# reducible polynomials will be ignored.
db_object = db.DeBruijnPoly('1011011')

# .polys property provides read-only access to the polynomials it was seeded with.
print 'created object, seeded with polynomials:'
print '  {}'.format(db_object.polys)

# .degree property provides easy access to the order of de bruijn sequences generated.
print 'object will generate debruijn sequences of order {}'.format(db_object.degree)

# .initial_state property allows de bruijn sequences to start on a specific state.
# defaults to an all-zero state.
print 'current initial state:'
print'  {}'.format(db_object.initial_state)

# .log property will cause object to time itself in each phase, as well as print it to console
# defaults to False
print 'object is logging execution time: {}'.format(db_object.log)

# .single_mode property will cause object to stop once a single de bruijn sequence can be
# generated. might or might not cut execution time drastically
print 'object is generating only one sequence: {}'.format(db_object.single_mode)

# object must be initialized first with .initialize() method. .initialized property indicates
# whether initialization process has finished or not
print 'object is initialized: {}'.format(db_object.initialized)
print '...'
db_object.initialize()
print 'object is initialized: {}'.format(db_object.initialized)

# once object is initialized, its adjacency matrix may be accessed to calculate how many
# sequences may be generated.
print 'number of db seqs: {}'.format(db_object.adjacency_matrix[1:, 1:].det())

# .ready property will indicate if object is currently ready to generate de bruijn sequences
print 'object is ready to generate sequence(s): {}'.format(db_object.ready)

# .anf property describes the algebraic normal form of the current sequence's modified feedback function
print 'algebraic normal form: {}'.format(db_object.anf)
print 'degree of ANF: {}'.format(db_object.anf_degree)

# .bits() method exposes a generator mechanism that outputs bits of the sequence one-by-one
# print 'first 10 bits of sequence:'
# ctr = 10
# for bit in db_object.bits():
#     print bit,
#     ctr -= 1
#     if ctr == 0:
#         break
# print

# fast-forward and exhaust the sequence
# for bit in db_object.bits():
#     pass

# once sequence is exhausted, object must be signalled to move on to the next de bruijn sequence
# by using .refresh_bits() method. this method passes -1 once all sequences are exhausted.
print 'object is ready to generate sequence(s): {}'.format(db_object.ready)
# print '.refresh_bits() returns: {}'.format(db_object.refresh_bits())

# .sequences() method exposes a generator mechanism that wraps around the .bits() generator
# and automates the .refresh_bits() process
print 'next sequence:'
for seq in db_object.sequences():
    print 'algebraic normal form: {}'.format(db_object.anf)
    print '  {}'.format(seq)
    break
