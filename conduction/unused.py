######## This file contains code written then removed from the main package. It may be useful in the future.
######## TB 10/27/16

## parallel tube generation (creation.Grid2D_onlat.init)

# if parallel:
#     comm = MPI.COMM_WORLD
#     logging.info("Setting up grid and tubes in parallel")
#     self.size = grid_size
#     if tube_length > grid_size:
#         logging.error('Nanotube is too large for grid')
#         raise SystemExit
#     # self holds all positions, on core 0 only
#     if rank == 0:
#         self.tube_coords = []
#         self.tube_coords_l = []
#         self.tube_coords_r = []
#         self.tube_centers = []
#         self.theta = []
#     else:
#         self.tube_coords = None
#         self.tube_coords_l = None
#         self.tube_coords_r = None
#         self.tube_centers = None
#         self.theta = None
#     self.tube_radius = tube_radius
#     counter = 0  # counts num of non-unique tubes replaced
#     whole_iterations = num_tubes / size
#     partial_iteration_num = num_tubes % size
#     if tube_radius == 0:
#         logging.info("Zero tube radius given. Tubes will have no volume.")
#         fill_fract = tube_length * float(num_tubes) / grid_size ** 2
#         logging.info("Filling fraction is %.2f %%" % (fill_fract * 100.0))
#         save_fill_frac(plot_save_dir, fill_fract)
#         comm.Barrier()
#         if num_tubes > 0:  # tubes exist
#             # generate all tubes on all cores, no checking yet
#             # whole iterations
#             for i in range(whole_iterations):
#                 x_l, y_l, x_r, y_r, x_c, y_c, theta = self.generate_2d_tube(tube_length, orientation,
#                                                                             tube_radius)
#                 # print('core %d here' % rank)
#                 local_tube_centers = [x_c, y_c]
#                 local_tube_coords = [x_l, y_l, x_r, y_r]
#                 local_tube_coords_l = [x_l, y_l]
#                 local_tube_coords_r = [x_r, y_r]
#                 local_theta = theta
#                 logging.info('Generating tube %d in parallel...' % (i * size))
#                 comm.Barrier()
#                 if rank == 0:
#                     # add work that was already done on core 0 to main lists
#                     self.tube_centers.append(local_tube_centers)
#                     self.tube_coords.append(local_tube_coords)
#                     self.tube_coords_l.append(local_tube_coords_l)
#                     self.tube_coords_r.append(local_tube_coords_r)
#                     self.theta.append(local_theta)
#                     for z in range(1, size):
#                         local_tube_centers_temp = comm.recv(source=z, tag=1)
#                         self.tube_centers.append(local_tube_centers_temp)
#                         local_tube_coords_temp = comm.recv(source=z, tag=2)
#                         # print local_tube_coords_temp
#                         self.tube_coords.append(local_tube_coords_temp)
#                         local_tube_coords_l_temp = comm.recv(source=z, tag=3)
#                         self.tube_coords_l.append(local_tube_coords_l_temp)
#                         local_tube_coords_r_temp = comm.recv(source=z, tag=4)
#                         self.tube_coords_r.append(local_tube_coords_r_temp)
#                         local_theta_temp = comm.recv(source=z, tag=5)
#                         self.theta.append(local_theta_temp)
#                 else:
#                     comm.send(local_tube_centers, dest=0, tag=1)
#                     comm.send(local_tube_coords, dest=0, tag=2)
#                     comm.send(local_tube_coords_l, dest=0, tag=3)
#                     comm.send(local_tube_coords_r, dest=0, tag=4)
#                     comm.send(local_theta, dest=0, tag=5)
#                 comm.Barrier()
#
#                 # code above checked and works for sending!
#                 # now send master *self* lists to all cores for checks
#                 # *self* will be stored on all cores (identical) and will hold the final values
#                 # A BCAST MUST SEND TO THE SAME VARIABLE ON ALL CORES, NOT A DIFFERENT ONE!
#                 self.tube_coords = comm.bcast(self.tube_coords, root=0)
#                 # print 'first bcast'
#                 # print self.tube_coords, rank
#                 uni_flag = self.check_tube_unique(self.tube_coords, True, rank,
#                                                   size)
#                 # print uni_flag, rank
#                 if rank != 0:
#                     comm.send(uni_flag, dest=0, tag=rank)
#                     uni_list = None
#                 else:
#                     uni_list = []  # stores list of non-unique tubes
#                     if not uni_flag:  # check on core 0
#                         uni_list.append(0)
#                     for z in range(1, size):
#                         uni_flag_temp = comm.recv(source=z, tag=z)
#                         if not uni_flag_temp:  # tube is not unique
#                             uni_list.append(z)
#                             # print uni_list
#                 uni_list = comm.bcast(uni_list, root=0)
#
#                 comm.Barrier()
#                 # uni_list holds which cores need to generate new tubes
#                 while uni_list:  # when it's not empty
#                     if rank in uni_list:
#                         counter += 1
#                         # generate new tube
#                         x_l, y_l, x_r, y_r, x_c, y_c, theta = self.generate_2d_tube(tube_length, orientation,
#                                                                                     tube_radius)
#                         # print('core %d here' % rank)
#                         local_tube_centers = [x_c, y_c]
#                         local_tube_coords = [x_l, y_l, x_r, y_r]
#                         local_tube_coords_l = [x_l, y_l]
#                         local_tube_coords_r = [x_r, y_r]
#                         local_theta = theta
#                     # comm.Barrier()
#                     # send to core 0
#                     if rank == 0:
#                         if rank in uni_list:
#                             temp = size - rank  # for rank 0 only
#                             # print temp
#                             self.tube_centers[-temp] = local_tube_centers
#                             self.tube_coords[-temp] = local_tube_coords
#                             self.tube_coords_l[-temp] = local_tube_coords_l
#                             self.tube_coords_r[-temp] = local_tube_coords_r
#                             self.theta[-temp] = local_theta
#                         for z in range(1, size):
#                             if z in uni_list:
#                                 temp = size - z
#                                 local_tube_centers_temp = comm.recv(source=z, tag=1 * z)
#                                 self.tube_centers[-temp] = local_tube_centers_temp
#                                 local_tube_coords_temp = comm.recv(source=z, tag=2 * z)
#                                 # print local_tube_coords_temp
#                                 self.tube_coords[-temp] = local_tube_coords_temp
#                                 local_tube_coords_l_temp = comm.recv(source=z, tag=3 * z)
#                                 self.tube_coords_l[-temp] = local_tube_coords_l_temp
#                                 local_tube_coords_r_temp = comm.recv(source=z, tag=4 * z)
#                                 self.tube_coords_r[-temp] = local_tube_coords_r_temp
#                                 local_theta_temp = comm.recv(source=z, tag=5 * z)
#                                 self.theta[-temp] = local_theta_temp
#                     else:
#                         if rank in uni_list:
#                             comm.send(local_tube_centers, dest=0, tag=1 * rank)
#                             comm.send(local_tube_coords, dest=0, tag=2 * rank)
#                             comm.send(local_tube_coords_l, dest=0, tag=3 * rank)
#                             comm.send(local_tube_coords_r, dest=0, tag=4 * rank)
#                             comm.send(local_theta, dest=0, tag=5 * rank)
#                     comm.Barrier()
#                     # regenerate uni_list
#                     self.tube_coords = comm.bcast(self.tube_coords, root=0)
#                     uni_flag = self.check_tube_unique(self.tube_coords, True, rank, size)
#                     if rank != 0:
#                         comm.send(uni_flag, dest=0, tag=rank)
#                         uni_list = None
#                     else:
#                         uni_list = []  # stores list of non-unique tubes
#                         if not uni_flag:  # check on core 0
#                             uni_list.append(0)
#                         for z in range(1, size):
#                             uni_flag_temp = comm.recv(source=z, tag=z)
#                             if not uni_flag_temp:  # tube is not unique
#                                 uni_list.append(z)
#                     uni_list = comm.bcast(uni_list, root=0)
#                     # end of regenerate tube while loop
#                     comm.Barrier()
#             # now handle any remainding tubes on core 0
#             if rank == 0:
#                 for c in range(partial_iteration_num):
#                     x_l, y_l, x_r, y_r, x_c, y_c, theta = self.generate_2d_tube(tube_length, orientation,
#                                                                                 tube_radius)
#                     self.tube_centers.append([x_c, y_c])
#                     self.tube_coords.append([x_l, y_l, x_r, y_r])
#                     self.tube_coords_l.append([x_l, y_l])
#                     self.tube_coords_r.append([x_r, y_r])
#                     self.theta.append(theta)
#                     logging.info('Generating tube %d on core 0...' % (i * size + c + 1))
#                     uni_flag = self.check_tube_unique(
#                         self.tube_coords, False)  # ensures no endpoints, left or right, are in the same spot
#                     while not uni_flag:
#                         counter += 1
#                         self.tube_centers.pop()
#                         self.tube_coords.pop()
#                         self.tube_coords_l.pop()
#                         self.tube_coords_r.pop()
#                         self.theta.pop()
#                         x_l, y_l, x_r, y_r, x_c, y_c, theta = self.generate_2d_tube(tube_length, orientation,
#                                                                                     tube_radius)
#                         self.tube_centers.append([x_c, y_c])
#                         self.tube_coords.append([x_l, y_l, x_r, y_r])
#                         self.tube_coords_l.append([x_l, y_l])
#                         self.tube_coords_r.append([x_r, y_r])
#                         self.theta.append(theta)
#                         uni_flag = self.check_tube_unique(self.tube_coords, False)
#             comm.Barrier()
#             # now that all tubes are unique, broadcast data to all cores
#             logging.info("Corrected %d overlapping tube endpoints" % counter)
#             self.tube_centers = comm.bcast(self.tube_centers, root=0)
#             self.tube_coords = comm.bcast(self.tube_coords, root=0)
#             self.tube_coords_l = comm.bcast(self.tube_coords_l, root=0)
#             self.tube_coords_r = comm.bcast(self.tube_coords_r, root=0)
#             self.theta = comm.bcast(self.theta, root=0)
#             comm.Barrier()
#             self.tube_check_l, self.tube_check_r, self.tube_check_bd = self.generate_tube_check_array_2d()
#     else:
#         logging.info("Non-zero tube radius given. Tubes will have excluded volume.")
#         if rank == 0:
#             self.tube_squares = []
#         else:
#             self.tube_squares = None
#         comm.Barrier()
#         if num_tubes > 0:  # tubes exist
#             # generate all tubes on all cores, no checking yet
#             # whole iterations
#             l_d = tube_length / (2 * tube_radius)
#             logging.info("L/D is %.4f." % l_d)
#             if num_tubes > 0:  # tubes exist
#                 for i in range(whole_iterations):  # currently no mean dist used, ADD LATER?
#                     x_l, y_l, x_r, y_r, x_c, y_c, theta = self.generate_2d_tube(tube_length, orientation,
#                                                                                 tube_radius)
#                     tube_squares = self.find_squares([x_l, y_l], [x_r, y_r], tube_radius)
#                     local_tube_centers = [x_c, y_c]
#                     local_tube_coords = [x_l, y_l, x_r, y_r]
#                     local_tube_coords_l = [x_l, y_l]
#                     local_tube_coords_r = [x_r, y_r]
#                     local_tube_squares = tube_squares
#                     local_theta = theta
#                     logging.info('Generating tube %d in parallel...' % (i * size))
#                     comm.Barrier()
#                     if rank == 0:
#                         # add work that was already done on core 0 to main lists
#                         self.tube_centers.append(local_tube_centers)
#                         self.tube_coords.append(local_tube_coords)
#                         self.tube_coords_l.append(local_tube_coords_l)
#                         self.tube_coords_r.append(local_tube_coords_r)
#                         self.theta.append(local_theta)
#                         self.tube_squares.append(local_tube_squares)
#                         for z in range(1, size):
#                             local_tube_centers_temp = comm.recv(source=z, tag=1)
#                             self.tube_centers.append(local_tube_centers_temp)
#                             local_tube_coords_temp = comm.recv(source=z, tag=2)
#                             self.tube_coords.append(local_tube_coords_temp)
#                             local_tube_coords_l_temp = comm.recv(source=z, tag=3)
#                             self.tube_coords_l.append(local_tube_coords_l_temp)
#                             local_tube_coords_r_temp = comm.recv(source=z, tag=4)
#                             self.tube_coords_r.append(local_tube_coords_r_temp)
#                             local_theta_temp = comm.recv(source=z, tag=5)
#                             self.theta.append(local_theta_temp)
#                             local_tube_squares_temp = comm.recv(source=z, tag=6)
#                             self.tube_squares.append(local_tube_squares_temp)
#                     else:
#                         comm.send(local_tube_centers, dest=0, tag=1)
#                         comm.send(local_tube_coords, dest=0, tag=2)
#                         comm.send(local_tube_coords_l, dest=0, tag=3)
#                         comm.send(local_tube_coords_r, dest=0, tag=4)
#                         comm.send(local_theta, dest=0, tag=5)
#                         comm.send(local_tube_squares, dest=0, tag=6)
#                     comm.Barrier()
#                     # code above checked and works for sending!
#                     # now send master *self* lists to all cores for checks
#                     # *self* will be stored on core 0 and will hold the final values
#                     # order starts to matter at this point!
#                     # A BCAST MUST SEND TO THE SAME VARIABLE ON ALL CORES, NOT A DIFFERENT ONE!
#                     self.tube_squares = comm.bcast(self.tube_squares, root=0)
#                     comm.Barrier()
#                     uni_flag = self.check_tube_and_vol_unique(self.tube_squares, True, rank,
#                                                               size)  # checks ends and volume
#                     comm.Barrier()
#
#                     if rank != 0:
#                         comm.send(uni_flag, dest=0, tag=rank)
#                         uni_list = None
#                     else:
#                         uni_list = []  # stores list of non-unique tubes
#                         if not uni_flag:  # check on core 0
#                             uni_list.append(0)
#                         for z in range(1, size):
#                             uni_flag_temp = comm.recv(source=z, tag=z)
#                             if not uni_flag_temp:  # tube is not unique
#                                 uni_list.append(z)
#
#                     comm.Barrier()
#                     uni_list = comm.bcast(uni_list, root=0)
#                     comm.Barrier()
#                     # uni_list holds which cores need to generate new tubes
#                     while uni_list:
#                         if rank in uni_list:
#                             counter += 1
#                             # generate new tube
#                             x_l, y_l, x_r, y_r, x_c, y_c, theta = self.generate_2d_tube(tube_length,
#                                                                                         orientation,
#                                                                                         tube_radius)
#                             tube_squares = self.find_squares([x_l, y_l], [x_r, y_r], tube_radius)
#                             # print('core %d here' % rank)
#                             local_tube_centers = [x_c, y_c]
#                             local_tube_coords = [x_l, y_l, x_r, y_r]
#                             local_tube_coords_l = [x_l, y_l]
#                             local_tube_coords_r = [x_r, y_r]
#                             local_theta = theta
#                             local_tube_squares = tube_squares
#                         comm.Barrier()
#                         # send to core 0
#                         if rank == 0:
#                             if rank in uni_list:
#                                 temp = size - rank  # SINCE INDEXING IS REVERSED, for rank 0 only
#                                 self.tube_centers[-temp] = local_tube_centers
#                                 self.tube_coords[-temp] = local_tube_coords
#                                 self.tube_coords_l[-temp] = local_tube_coords_l
#                                 self.tube_coords_r[-temp] = local_tube_coords_r
#                                 self.theta[-temp] = local_theta
#                                 self.tube_squares[-temp] = local_tube_squares
#                             for z in range(1, size):
#                                 if z in uni_list:
#                                     temp = size - z
#                                     local_tube_centers_temp = comm.recv(source=z, tag=1)
#                                     self.tube_centers[-temp] = local_tube_centers_temp
#                                     local_tube_coords_temp = comm.recv(source=z, tag=2)
#                                     self.tube_coords[-temp] = local_tube_coords_temp
#                                     local_tube_coords_l_temp = comm.recv(source=z, tag=3)
#                                     self.tube_coords_l[-temp] = local_tube_coords_l_temp
#                                     local_tube_coords_r_temp = comm.recv(source=z, tag=4)
#                                     self.tube_coords_r[-temp] = local_tube_coords_r_temp
#                                     local_theta_temp = comm.recv(source=z, tag=5)
#                                     self.theta[-temp] = local_theta_temp
#                                     local_tube_squares_temp = comm.recv(source=z, tag=6)
#                                     self.tube_squares[-temp] = local_tube_squares_temp
#                                     # print 'Received for item %d on core %d' % (z, rank)
#                         else:
#                             if rank in uni_list:
#                                 comm.send(local_tube_centers, dest=0, tag=1)
#                                 comm.send(local_tube_coords, dest=0, tag=2)
#                                 comm.send(local_tube_coords_l, dest=0, tag=3)
#                                 comm.send(local_tube_coords_r, dest=0, tag=4)
#                                 comm.send(local_theta, dest=0, tag=5)
#                                 comm.send(local_tube_squares, dest=0, tag=6)
#                                 # print 'Sent from core %d' % rank
#                         # regenerate uni_list
#                         comm.Barrier()
#                         self.tube_squares = comm.bcast(self.tube_squares, root=0)  # sends to all but core 0
#                         comm.Barrier()
#                         uni_flag = self.check_tube_and_vol_unique(self.tube_squares, True, rank,
#                                                                   size)  # checks ends and volume
#                         comm.Barrier()
#                         if rank != 0:
#                             comm.send(uni_flag, dest=0, tag=rank)
#                             uni_list = None
#                         else:
#                             # print uni_flag
#                             uni_list = []  # stores list of non-unique tubes
#                             if not uni_flag:  # check on core 0
#                                 uni_list.append(0)
#                             for z in range(1, size):
#                                 uni_flag_temp = comm.recv(source=z, tag=z)
#                                 if not uni_flag_temp:  # tube is not unique
#                                     uni_list.append(z)
#                         uni_list = comm.bcast(uni_list, root=0)
#                         # end of regenerate tube while loop
#                         comm.Barrier()
#                 # now handle any remaining tubes on core 0
#                 if rank == 0:
#                     for c in range(partial_iteration_num):
#                         x_l, y_l, x_r, y_r, x_c, y_c, theta = self.generate_2d_tube(tube_length, orientation,
#                                                                                     tube_radius)
#                         tube_squares = self.find_squares([x_l, y_l], [x_r, y_r], tube_radius)
#                         self.tube_centers.append([x_c, y_c])
#                         self.tube_coords.append([x_l, y_l, x_r, y_r])
#                         self.tube_coords_l.append([x_l, y_l])
#                         self.tube_coords_r.append([x_r, y_r])
#                         self.tube_squares.append(tube_squares)
#                         self.theta.append(theta)
#                         logging.info('Generating tube %d on core 0...' % (i * size + c + 1))
#                         uni_flag = self.check_tube_and_vol_unique(self.tube_squares, False)
#                         while not uni_flag:
#                             counter += 1
#                             self.tube_centers.pop()
#                             self.tube_coords.pop()
#                             self.tube_coords_l.pop()
#                             self.tube_coords_r.pop()
#                             self.tube_squares.pop()
#                             self.theta.pop()
#                             x_l, y_l, x_r, y_r, x_c, y_c, theta = self.generate_2d_tube(tube_length,
#                                                                                         orientation,
#                                                                                         tube_radius)
#                             tube_squares = self.find_squares([x_l, y_l], [x_r, y_r], tube_radius)
#                             self.theta.append(theta)
#                             self.tube_centers.append([x_c, y_c])
#                             self.tube_coords.append([x_l, y_l, x_r, y_r])
#                             self.tube_coords_l.append([x_l, y_l])
#                             self.tube_coords_r.append([x_r, y_r])
#                             self.tube_squares.append(tube_squares)
#                             uni_flag = self.check_tube_and_vol_unique(self.tube_squares, False)
#                 comm.Barrier()
#                 # now that all tubes are unique, broadcast data to all cores
#                 logging.info("Corrected %d overlapping tube endpoints" % counter)
#                 self.tube_centers = comm.bcast(self.tube_centers, root=0)
#                 self.tube_coords = comm.bcast(self.tube_coords, root=0)
#                 self.tube_coords_l = comm.bcast(self.tube_coords_l, root=0)
#                 self.tube_coords_r = comm.bcast(self.tube_coords_r, root=0)
#                 self.theta = comm.bcast(self.theta, root=0)
#                 self.tube_squares = comm.bcast(self.tube_squares, root=0)
#                 comm.Barrier()
#                 # this stuff can be done on all cores since they have all the data now,
#                 # cause why the hell not?
#                 # get number of squares filled
#                 cube_count = 0  # each cube has area 1
#                 for i in range(len(self.tube_squares)):
#                     cube_count += len(self.tube_squares[i])
#                 fill_fract = float(cube_count) * 2.0 * tube_radius / grid_size ** 2
#                 # each cube has area 1, times the tube radius (important if not 1)
#                 logging.info("Filling fraction is %.2f %%" % (fill_fract * 100.0))
#                 save_fill_frac(plot_save_dir, fill_fract)
#                 self.tube_check_l, self.tube_check_r, self.tube_check_bd = self.generate_tube_check_array_2d()
#                 self.tube_check_bd_vol, self.tube_check_index = self.generate_vol_check_array_2d()

## parallel tube generation (creation.Grid3D_onlat.init)

# if parallel:
#     comm = MPI.COMM_WORLD
#     logging.info("Setting up grid and tubes in parallel")
#     if tube_length > grid_size:
#         logging.error('Nanotube is too large for grid')
#         raise SystemExit
#     # self holds all positions, on core 0 only
#     if rank == 0:
#         self.tube_coords = []
#         self.tube_coords_l = []
#         self.tube_coords_r = []
#         self.tube_centers = []
#         self.theta = []
#         self.phi = []
#     else:
#         self.tube_coords = None
#         self.tube_coords_l = None
#         self.tube_coords_r = None
#         self.tube_centers = None
#         self.theta = None
#         self.phi = None
#     counter = 0  # counts num of non-unique tubes replaced
#     whole_iterations = num_tubes / size
#     partial_iteration_num = num_tubes % size
#     if tube_radius == 0:
#         logging.info("Zero tube radius given. Tubes will have no volume.")
#         fill_fract = tube_length * float(num_tubes) / grid_size ** 3
#         logging.info("Filling fraction is %.2f %%" % (fill_fract * 100.0))
#         save_fill_frac(plot_save_dir, fill_fract)
#         comm.Barrier()
#         if num_tubes > 0:  # tubes exist
#             # generate all tubes on all cores, no checking yet
#             # whole iterations
#             for i in range(whole_iterations):
#                 x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c, theta, phi = self.generate_3d_tube(tube_length,
#                                                                                                 orientation,
#                                                                                                 tube_radius)
#                 local_tube_centers = [x_c, y_c, z_c]
#                 local_tube_coords = [x_l, y_l, z_l, x_r, y_r, z_r]
#                 local_tube_coords_l = [x_l, y_l, z_l]
#                 local_tube_coords_r = [x_r, y_r, z_r]
#                 local_theta = theta
#                 local_phi = phi
#                 logging.info('Generating tube %d in parallel...' % (i * size))
#                 comm.Barrier()
#                 if rank == 0:
#                     # add work that was already done on core 0 to main lists
#                     self.tube_centers.append(local_tube_centers)
#                     self.tube_coords.append(local_tube_coords)
#                     self.tube_coords_l.append(local_tube_coords_l)
#                     self.tube_coords_r.append(local_tube_coords_r)
#                     self.theta.append(local_theta)
#                     self.phi.append(local_phi)
#                     for z in range(1, size):
#                         local_tube_centers_temp = comm.recv(source=z, tag=1)
#                         self.tube_centers.append(local_tube_centers_temp)
#                         local_tube_coords_temp = comm.recv(source=z, tag=2)
#                         self.tube_coords.append(local_tube_coords_temp)
#                         local_tube_coords_l_temp = comm.recv(source=z, tag=3)
#                         self.tube_coords_l.append(local_tube_coords_l_temp)
#                         local_tube_coords_r_temp = comm.recv(source=z, tag=4)
#                         self.tube_coords_r.append(local_tube_coords_r_temp)
#                         local_theta_temp = comm.recv(source=z, tag=5)
#                         self.theta.append(local_theta_temp)
#                         local_phi_temp = comm.recv(source=z, tag=6)
#                         self.phi.append(local_phi_temp)
#                 else:
#                     comm.send(local_tube_centers, dest=0, tag=1)
#                     comm.send(local_tube_coords, dest=0, tag=2)
#                     comm.send(local_tube_coords_l, dest=0, tag=3)
#                     comm.send(local_tube_coords_r, dest=0, tag=4)
#                     comm.send(local_theta, dest=0, tag=5)
#                     comm.send(local_phi, dest=0, tag=6)
#                 comm.Barrier()
#                 # code above checked and works for sending!
#                 # now send master *self* lists to all cores for checks
#                 # *self* will be stored on all cores (identical) and will hold the final values
#                 # A BCAST MUST SEND TO THE SAME VARIABLE ON ALL CORES, NOT A DIFFERENT ONE!
#                 self.tube_coords = comm.bcast(self.tube_coords, root=0)
#                 uni_flag = self.check_tube_unique(self.tube_coords, True, rank, size)
#                 # print uni_flag, rank
#                 if rank != 0:
#                     comm.send(uni_flag, dest=0, tag=rank)
#                     uni_list = None
#                 else:
#                     uni_list = []  # stores list of non-unique tubes
#                     if not uni_flag:  # check on core 0
#                         uni_list.append(0)
#                     for z in range(1, size):
#                         uni_flag_temp = comm.recv(source=z, tag=z)
#                         if not uni_flag_temp:  # tube is not unique
#                             uni_list.append(z)
#                             # print uni_list
#                 uni_list = comm.bcast(uni_list, root=0)
#                 # print uni_list
#                 comm.Barrier()
#                 # print 'here' # good
#                 # uni_list holds which cores need to generate new tubes
#                 while uni_list:  # when it's not empty
#                     if rank in uni_list:
#                         counter += 1
#                         # generate new tube
#                         x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c, theta, phi = self.generate_3d_tube(
#                             tube_length,
#                             orientation,
#                             tube_radius)
#                         local_tube_centers = [x_c, y_c, z_c]
#                         local_tube_coords = [x_l, y_l, z_l, x_r, y_r, z_r]
#                         local_tube_coords_l = [x_l, y_l, z_l]
#                         local_tube_coords_r = [x_r, y_r, z_r]
#                         local_theta = theta
#                         local_phi = phi
#                     comm.Barrier()
#                     # send to core 0
#                     if rank == 0:
#                         if rank in uni_list:
#                             temp = size - rank  # SINCE INDEXING IS REVERSED, for rank 0 only
#                             # print temp
#                             self.tube_centers[-temp] = local_tube_centers
#                             self.tube_coords[-temp] = local_tube_coords
#                             self.tube_coords_l[-temp] = local_tube_coords_l
#                             self.tube_coords_r[-temp] = local_tube_coords_r
#                             self.theta[-temp] = local_theta
#                             self.phi[-temp] = local_phi
#                         for z in range(1, size):
#                             if z in uni_list:
#                                 temp = size - z
#                                 local_tube_centers_temp = comm.recv(source=z, tag=1 * z)
#                                 self.tube_centers[-temp] = local_tube_centers_temp
#                                 local_tube_coords_temp = comm.recv(source=z, tag=2 * z)
#                                 # print local_tube_coords_temp
#                                 self.tube_coords[-temp] = local_tube_coords_temp
#                                 local_tube_coords_l_temp = comm.recv(source=z, tag=3 * z)
#                                 self.tube_coords_l[-temp] = local_tube_coords_l_temp
#                                 local_tube_coords_r_temp = comm.recv(source=z, tag=4 * z)
#                                 self.tube_coords_r[-temp] = local_tube_coords_r_temp
#                                 local_theta_temp = comm.recv(source=z, tag=5 * z)
#                                 self.theta[-temp] = local_theta_temp
#                                 local_phi_temp = comm.recv(source=z, tag=6 * z)
#                                 self.phi[-temp] = local_phi_temp
#                     else:
#                         if rank in uni_list:
#                             comm.send(local_tube_centers, dest=0, tag=1 * rank)
#                             comm.send(local_tube_coords, dest=0, tag=2 * rank)
#                             comm.send(local_tube_coords_l, dest=0, tag=3 * rank)
#                             comm.send(local_tube_coords_r, dest=0, tag=4 * rank)
#                             comm.send(local_theta, dest=0, tag=5 * rank)
#                             comm.send(local_phi, dest=0, tag=6 * rank)
#                     comm.Barrier()
#                     # regenerate uni_list
#                     # print 'here' # good
#                     # if rank == 0:
#                     # print checker_tube_coords
#                     self.tube_coords = comm.bcast(self.tube_coords, root=0)
#                     uni_flag = self.check_tube_unique(self.tube_coords, True, rank, size)
#
#                     if rank != 0:
#                         comm.send(uni_flag, dest=0, tag=rank)
#                         uni_list = None
#                     else:
#                         # print uni_flag
#                         uni_list = []  # stores list of non-unique tubes
#                         if not uni_flag:  # check on core 0
#                             uni_list.append(0)
#                         for z in range(1, size):
#                             uni_flag_temp = comm.recv(source=z, tag=z)
#                             if not uni_flag_temp:  # tube is not unique
#                                 uni_list.append(z)
#                     uni_list = comm.bcast(uni_list, root=0)
#                     comm.Barrier()
#                     # now handle any remainding tubes on core 0
#             if rank == 0:
#                 for c in range(partial_iteration_num):
#                     x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c, theta, phi = self.generate_3d_tube(tube_length,
#                                                                                                     orientation,
#                                                                                                     tube_radius)
#                     self.tube_centers.append([x_c, y_c, z_c])
#                     self.tube_coords.append([x_l, y_l, z_l, x_r, y_r, z_r])
#                     self.tube_coords_l.append([x_l, y_l, z_l])
#                     self.tube_coords_r.append([x_r, y_r, z_r])
#                     self.theta.append(theta)
#                     self.phi.append(phi)
#                     logging.info('Generating tube %d on core 0...' % (i * size + c + 1))
#                     uni_flag = self.check_tube_unique(
#                         self.tube_coords, False)  # ensures no endpoints, left or right, are in the same spot
#                     while not uni_flag:
#                         counter += 1
#                         self.tube_centers.pop()
#                         self.tube_coords.pop()
#                         self.tube_coords_l.pop()
#                         self.tube_coords_r.pop()
#                         self.theta.pop()
#                         self.phi.pop()
#                         x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c, theta, phi = self.generate_3d_tube(
#                             tube_length,
#                             orientation,
#                             tube_radius)
#                         self.tube_centers.append([x_c, y_c, z_c])
#                         self.tube_coords.append([x_l, y_l, z_l, x_r, y_r, z_r])
#                         self.tube_coords_l.append([x_l, y_l, z_l])
#                         self.tube_coords_r.append([x_r, y_r, z_r])
#                         self.theta.append(theta)
#                         self.phi.append(phi)
#                         uni_flag = self.check_tube_unique(self.tube_coords, False)
#             comm.Barrier()
#             # now that all tubes are unique, broadcast data to all cores
#             self.tube_centers = comm.bcast(self.tube_centers, root=0)
#             self.tube_coords = comm.bcast(self.tube_coords, root=0)
#             self.tube_coords_l = comm.bcast(self.tube_coords_l, root=0)
#             self.tube_coords_r = comm.bcast(self.tube_coords_r, root=0)
#             self.theta = comm.bcast(self.theta, root=0)
#             self.phi = comm.bcast(self.phi, root=0)
#             comm.Barrier()
#         logging.info("Corrected %d overlapping tube endpoints" % counter)
#         logging.info("Tube generation complete")
#         self.tube_check_l, self.tube_check_r, self.tube_check_bd = self.generate_tube_check_array_3d()
#     else:
#         if rank == 0:
#             self.tube_squares = []
#         else:
#             self.tube_squares = None
#         comm.Barrier()
#         if num_tubes > 0:  # tubes exist
#             logging.info("Non-zero tube radius given. Tubes will have excluded volume.")
#             l_d = tube_length / (2 * tube_radius)
#             logging.info("L/D is %.4f." % l_d)
#             for i in range(whole_iterations):
#                 x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c, theta, phi = self.generate_3d_tube(tube_length,
#                                                                                                 orientation,
#                                                                                                 tube_radius)
#                 tube_squares = self.find_cubes([x_l, y_l, z_l], [x_r, y_r, z_r])
#                 local_tube_centers = [x_c, y_c, z_c]
#                 local_tube_coords = [x_l, y_l, z_l, x_r, y_r, z_r]
#                 local_tube_coords_l = [x_l, y_l, z_l]
#                 local_tube_coords_r = [x_r, y_r, z_r]
#                 local_theta = theta
#                 local_phi = phi
#                 local_tube_squares = tube_squares
#                 logging.info('Generating tube %d in parallel...' % (i * size))
#                 comm.Barrier()
#                 if rank == 0:
#                     # add work that was already done on core 0 to main lists
#                     self.tube_centers.append(local_tube_centers)
#                     self.tube_coords.append(local_tube_coords)
#                     self.tube_coords_l.append(local_tube_coords_l)
#                     self.tube_coords_r.append(local_tube_coords_r)
#                     self.theta.append(local_theta)
#                     self.phi.append(local_phi)
#                     self.tube_squares.append(local_tube_squares)
#                     for z in range(1, size):
#                         local_tube_centers_temp = comm.recv(source=z, tag=1)
#                         self.tube_centers.append(local_tube_centers_temp)
#                         local_tube_coords_temp = comm.recv(source=z, tag=2)
#                         self.tube_coords.append(local_tube_coords_temp)
#                         local_tube_coords_l_temp = comm.recv(source=z, tag=3)
#                         self.tube_coords_l.append(local_tube_coords_l_temp)
#                         local_tube_coords_r_temp = comm.recv(source=z, tag=4)
#                         self.tube_coords_r.append(local_tube_coords_r_temp)
#                         local_theta_temp = comm.recv(source=z, tag=5)
#                         self.theta.append(local_theta_temp)
#                         local_phi_temp = comm.recv(source=z, tag=6)
#                         self.phi.append(local_phi_temp)
#                         local_tube_squares_temp = comm.recv(source=z, tag=7)
#                         self.tube_squares.append(local_tube_squares_temp)
#                 else:
#                     comm.send(local_tube_centers, dest=0, tag=1)
#                     comm.send(local_tube_coords, dest=0, tag=2)
#                     comm.send(local_tube_coords_l, dest=0, tag=3)
#                     comm.send(local_tube_coords_r, dest=0, tag=4)
#                     comm.send(local_theta, dest=0, tag=5)
#                     comm.send(local_phi, dest=0, tag=6)
#                     comm.send(local_tube_squares, dest=0, tag=7)
#                 comm.Barrier()
#                 # code above checked and works for sending!
#                 # now send master *self* lists to all cores for checks
#                 # *self* will be stored on all cores (identical) and will hold the final values
#                 # A BCAST MUST SEND TO THE SAME VARIABLE ON ALL CORES, NOT A DIFFERENT ONE!
#                 self.tube_squares = comm.bcast(self.tube_squares, root=0)
#                 uni_flag = self.check_tube_and_vol_unique(self.tube_squares, True, rank, size)
#                 # print uni_flag, rank
#                 if rank != 0:
#                     comm.send(uni_flag, dest=0, tag=rank)
#                     uni_list = None
#                 else:
#                     uni_list = []  # stores list of non-unique tubes
#                     if not uni_flag:  # check on core 0
#                         uni_list.append(0)
#                     for z in range(1, size):
#                         uni_flag_temp = comm.recv(source=z, tag=z)
#                         if not uni_flag_temp:  # tube is not unique
#                             uni_list.append(z)
#                             # print uni_list
#                 uni_list = comm.bcast(uni_list, root=0)
#                 comm.Barrier()
#                 # uni_list holds which cores need to generate new tubes
#                 while uni_list:  # when it's not empty
#                     if rank in uni_list:
#                         counter += 1
#                         # generate new tube
#                         x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c, theta, phi = self.generate_3d_tube(
#                             tube_length,
#                             orientation,
#                             tube_radius)
#                         tube_squares = self.find_cubes([x_l, y_l, z_l], [x_r, y_r, z_r])
#                         local_tube_centers = [x_c, y_c, z_c]
#                         local_tube_coords = [x_l, y_l, z_l, x_r, y_r, z_r]
#                         local_tube_coords_l = [x_l, y_l, z_l]
#                         local_tube_coords_r = [x_r, y_r, z_r]
#                         local_theta = theta
#                         local_phi = phi
#                         local_tube_squares = tube_squares
#                     comm.Barrier()
#                     # send to core 0
#                     if rank == 0:
#                         if rank in uni_list:
#                             temp = size - rank  # SINCE INDEXING IS REVERSED, for rank 0 only
#                             # print temp
#                             self.tube_centers[-temp] = local_tube_centers
#                             self.tube_coords[-temp] = local_tube_coords
#                             self.tube_coords_l[-temp] = local_tube_coords_l
#                             self.tube_coords_r[-temp] = local_tube_coords_r
#                             self.theta[-temp] = local_theta
#                             self.phi[-temp] = local_phi
#                             self.tube_squares[-temp] = local_tube_squares
#                         for z in range(1, size):
#                             if z in uni_list:
#                                 temp = size - z
#                                 local_tube_centers_temp = comm.recv(source=z, tag=1 * z)
#                                 self.tube_centers[-temp] = local_tube_centers_temp
#                                 local_tube_coords_temp = comm.recv(source=z, tag=2 * z)
#                                 # print local_tube_coords_temp
#                                 self.tube_coords[-temp] = local_tube_coords_temp
#                                 local_tube_coords_l_temp = comm.recv(source=z, tag=3 * z)
#                                 self.tube_coords_l[-temp] = local_tube_coords_l_temp
#                                 local_tube_coords_r_temp = comm.recv(source=z, tag=4 * z)
#                                 self.tube_coords_r[-temp] = local_tube_coords_r_temp
#                                 local_theta_temp = comm.recv(source=z, tag=5 * z)
#                                 self.theta[-temp] = local_theta_temp
#                                 local_phi_temp = comm.recv(source=z, tag=6 * z)
#                                 self.phi[-temp] = local_phi_temp
#                                 local_tube_squares_temp = comm.recv(source=z, tag=7 * z)
#                                 self.tube_squares[-temp] = local_tube_squares_temp
#                     else:
#                         if rank in uni_list:
#                             comm.send(local_tube_centers, dest=0, tag=1 * rank)
#                             comm.send(local_tube_coords, dest=0, tag=2 * rank)
#                             comm.send(local_tube_coords_l, dest=0, tag=3 * rank)
#                             comm.send(local_tube_coords_r, dest=0, tag=4 * rank)
#                             comm.send(local_theta, dest=0, tag=5 * rank)
#                             comm.send(local_phi, dest=0, tag=6 * rank)
#                             comm.send(local_tube_squares, dest=0, tag=7 * rank)
#                     comm.Barrier()
#                     # regenerate uni_list
#                     # print 'here' # good
#                     # if rank == 0:
#                     # print checker_tube_coords
#                     self.tube_squares = comm.bcast(self.tube_squares, root=0)
#                     uni_flag = self.check_tube_and_vol_unique(self.tube_squares, True, rank, size)
#                     # print uni_flag, rank
#                     if rank != 0:
#                         comm.send(uni_flag, dest=0, tag=rank)
#                         uni_list = None
#                     else:
#                         # print uni_flag
#                         uni_list = []  # stores list of non-unique tubes
#                         if not uni_flag:  # check on core 0
#                             uni_list.append(0)
#                         for z in range(1, size):
#                             uni_flag_temp = comm.recv(source=z, tag=z)
#                             if not uni_flag_temp:  # tube is not unique
#                                 uni_list.append(z)
#                     uni_list = comm.bcast(uni_list, root=0)
#                     comm.Barrier()
#                     # now handle any remainding tubes on core 0
#             if rank == 0:
#                 for c in range(partial_iteration_num):
#                     x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c, theta, phi = self.generate_3d_tube(tube_length,
#                                                                                                     orientation,
#                                                                                                     tube_radius)
#                     tube_squares = self.find_cubes([x_l, y_l, z_l], [x_r, y_r, z_r])
#                     self.tube_centers.append([x_c, y_c, z_c])
#                     self.tube_coords.append([x_l, y_l, z_l, x_r, y_r, z_r])
#                     self.tube_coords_l.append([x_l, y_l, z_l])
#                     self.tube_coords_r.append([x_r, y_r, z_r])
#                     self.theta.append(theta)
#                     self.phi.append(phi)
#                     self.tube_squares.append(tube_squares)
#                     logging.info('Generating tube %d on core 0...' % (i * size + c + 1))
#                     uni_flag = self.check_tube_and_vol_unique(self.tube_squares, False)
#                     while not uni_flag:
#                         counter += 1
#                         self.tube_centers.pop()
#                         self.tube_coords.pop()
#                         self.tube_coords_l.pop()
#                         self.tube_coords_r.pop()
#                         self.theta.pop()
#                         self.phi.pop()
#                         self.tube_squares.pop()
#                         x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c, theta, phi = self.generate_3d_tube(
#                             tube_length,
#                             orientation,
#                             tube_radius)
#                         tube_squares = self.find_cubes([x_l, y_l, z_l], [x_r, y_r, z_r])
#                         self.tube_centers.append([x_c, y_c, z_c])
#                         self.tube_coords.append([x_l, y_l, z_l, x_r, y_r, z_r])
#                         self.tube_coords_l.append([x_l, y_l, z_l])
#                         self.tube_coords_r.append([x_r, y_r, z_r])
#                         self.theta.append(theta)
#                         self.phi.append(phi)
#                         self.tube_squares.append(tube_squares)
#                         uni_flag = self.check_tube_and_vol_unique(self.tube_squares, False)
#             comm.Barrier()
#             # now that all tubes are unique, broadcast data to all cores
#             logging.info("Tube generation complete")
#             logging.info("Corrected %d overlapping tube endpoints" % counter)
#             self.tube_centers = comm.bcast(self.tube_centers, root=0)
#             self.tube_coords = comm.bcast(self.tube_coords, root=0)
#             self.tube_coords_l = comm.bcast(self.tube_coords_l, root=0)
#             self.tube_coords_r = comm.bcast(self.tube_coords_r, root=0)
#             self.theta = comm.bcast(self.theta, root=0)
#             self.phi = comm.bcast(self.phi, root=0)
#             self.tube_squares = comm.bcast(self.tube_squares, root=0)
#             comm.Barrier()
#         logging.info("Corrected %d overlapping tube endpoints" % counter)
#         # get number of squares filled
#         cube_count = 0  # each cube has volume 1
#         for i in range(len(self.tube_squares)):
#             cube_count += len(self.tube_squares[i])
#         fill_fract = float(cube_count) * 2.0 * tube_radius / grid_size ** 3
#         # each cube has area 1, times the tube radius (important if not 1)
#         logging.info("Filling fraction is %.2f %%" % (fill_fract * 100.0))
#         save_fill_frac(plot_save_dir, fill_fract)
#         self.tube_check_l, self.tube_check_r, self.tube_check_bd = self.generate_tube_check_array_3d()
#         self.tube_check_bd_vol, self.tube_check_index = self.generate_vol_check_array_3d(disable_func)
