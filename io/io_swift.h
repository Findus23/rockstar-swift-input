#ifndef _IO_SWIFT_H_
#define _IO_SWIFT_H_
#ifdef ENABLE_HDF5
#include <stdint.h>
#include "../particle.h"

void load_particles_swift(char *filename, struct particle **p, int64_t *num_p);

#endif /* ENABLE_HDF5 */
#endif /* _IO_SWIFT_H_ */
