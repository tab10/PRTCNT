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


def generate_bd_and_vol(self, tube_squares, theta, phi):
    """3D, This takes the bottom left squares from find_squares and creates volume and boundaries based
    on the tube radius, for one tube
    In 3D, there are top left, top right, bottom left, and bottom right boundaries to worry about too?
    Currently not used for anything
    """
    top = []
    bot = []
    left = []
    right = []
    l_end = tube_squares[0]
    r_end = tube_squares[-1]
    occupied_cubes = tube_squares  # this will be checked to ensure volume is exclusive too
    for i in range(1, self.tube_radius + 1):
        l_x_above = int(round(self.coord(i, theta + 90, phi + 90)[0] + l_end[0]))
        l_y_above = int(round(self.coord(i, theta + 90, phi + 90)[1] + l_end[1]))
        l_z_above = int(round(self.coord(i, theta + 90, phi + 90)[2] + l_end[2]))
        r_x_above = int(round(self.coord(i, theta + 90, phi + 90)[0] + r_end[0]))
        r_y_above = int(round(self.coord(i, theta + 90, phi + 90)[1] + r_end[1]))
        r_z_above = int(round(self.coord(i, theta + 90, phi + 90)[2] + r_end[2]))
        l_x_below = int(round(self.coord(i, theta - 90, phi - 90)[0] + l_end[0]))
        l_y_below = int(round(self.coord(i, theta - 90, phi - 90)[1] + l_end[1]))
        l_z_below = int(round(self.coord(i, theta - 90, phi - 90)[2] + l_end[2]))
        r_x_below = int(round(self.coord(i, theta - 90, phi - 90)[0] + r_end[0]))
        r_y_below = int(round(self.coord(i, theta - 90, phi - 90)[1] + r_end[1]))
        r_z_below = int(round(self.coord(i, theta - 90, phi - 90)[2] + r_end[2]))
        left.append([l_x_above, l_y_above, l_z_above])
        left.append([l_x_below, l_y_below, l_z_below])
        right.append([r_x_above, r_y_above, r_z_above])
        right.append([r_x_below, r_y_below, r_z_below])
    for i in range(len(tube_squares)):
        t_x = int(round(self.coord(self.tube_radius, theta + 90, phi + 90)[0] + tube_squares[i][0]))
        t_y = int(round(self.coord(self.tube_radius, theta + 90, phi + 90)[1] + tube_squares[i][1]))
        t_z = int(round(self.coord(self.tube_radius, theta + 90, phi + 90)[2] + tube_squares[i][2]))
        b_x = int(round(self.coord(self.tube_radius, theta - 90, phi - 90)[0] + tube_squares[i][0]))
        b_y = int(round(self.coord(self.tube_radius, theta - 90, phi - 90)[1] + tube_squares[i][1]))
        b_z = int(round(self.coord(self.tube_radius, theta - 90, phi - 90)[2] + tube_squares[i][2]))
        top.append([t_x, t_y, t_z])
        bot.append([b_x, b_y, b_z])
    total = top + bot + left + right
    return top, bot, left, right, total, occupied_cubes


## OLD BD COND 3D FUNCTION ##
# def apply_bd_cond_3d(grid, moves_3d, cur_pos, bound):
#     choices = generate_bd_choices_3d(grid, cur_pos, moves_3d, bound)
#     # pick random choice
#     new_pos = np.asarray(choices[np.random.randint(0, len(choices))])
#     final_pos = new_pos
# if rules_test:
#     choices = generate_periodic_bc_choices_3d(grid, cur_pos, moves_3d)
#     # pick random choice
#     new_pos = np.asarray(choices[np.random.randint(0, 6)])
#     final_pos = new_pos
# if cur_type == 10:  # reflective boundary
#     choices = generate_periodic_bc_choices_3d(grid, cur_pos, moves_3d)
#     # pick random choice
#     d_pos = np.asarray(moves_3d[np.random.randint(0, 6)])
#     candidate_pos = cur_pos + d_pos
#     if (candidate_pos[0] > grid.size) or (candidate_pos[0] < 0):
#         final_pos = cur_pos
#     else:
#         final_pos = candidate_pos
# elif cur_type == 20:  # periodic boundary
#     d_pos = np.asarray(moves_3d[np.random.randint(0, 6)])
#     candidate_pos = cur_pos + d_pos
#     if candidate_pos[1] > grid.size:
#         final_pos = [candidate_pos[0], 1, candidate_pos[2]]
#     elif candidate_pos[1] < 0:
#         final_pos = [candidate_pos[0], grid.size - 1, candidate_pos[2]]
#     elif candidate_pos[2] > grid.size:
#         final_pos = [candidate_pos[0], candidate_pos[1], 1]
#     elif candidate_pos[2] < 0:
#         final_pos = [candidate_pos[0], candidate_pos[1], grid.size - 1]
#     else:
#         final_pos = candidate_pos
# elif cur_type == 30:  # corner
#     d_pos = np.asarray(moves_3d[np.random.randint(0, 6)])
#     candidate_pos = cur_pos + d_pos
#     if (candidate_pos[0] > grid.size) or (candidate_pos[0] < 0) or (candidate_pos[1] > grid.size) or (candidate_pos[1] < 0) \
#             or (candidate_pos[2] > grid.size) or (candidate_pos[2] < 0):
#         final_pos = cur_pos
#     else:
#         final_pos = candidate_pos
# else:
#     kill()
# return final_pos


"""These were very inefficient and were replaced. TB 4/9/2017"""
# def generate_periodic_bc_choices_3d(grid, cur_pos, moves_3d):
#     """Returns a list of move choices for a walker on a periodic boundary"""
#     choices = cur_pos + moves_3d
#     # get edge locations
#     min_val = 0
#     max_val = grid.size  # not + 1 as in setup as we can walk on 0 or 100
#     # check choices for crossover
#     for i in range(len(choices)):
#         temp = choices[i]
#         for j in range(len(temp)):
#             if temp[j] < min_val:
#                 temp[j] = max_val
#             elif temp[j] > max_val:
#                 temp[j] = min_val
#         choices[i] = temp
#     return choices


# def generate_reflective_bc_choices_3d(grid, cur_pos, moves_3d):
#     """Returns a list of move choices for a walker on a reflective boundary"""
#     choices = cur_pos + moves_3d
#     # get edge locations
#     min_val = 0
#     max_val = grid.size  # not + 1 as in setup as we can walk on 0 or 100
#     # we are guaranteed that that coordinate outside bd is reflective because we're told that :P
#     for i in range(len(choices)):
#         temp = choices[i]
#         for j in range(len(temp)):
#             if (temp[j] < min_val) or (temp[j] > max_val):
#                 temp[j] = None
#                 break
#         if None in temp:
#             del choices[i]
#         else:
#             choices[i] = temp
#     return choices


# def generate_corner_bc_choices_3d(grid, cur_pos, moves_3d):
#     """Returns a list of move choices for a walker on a corner boundary"""
#     choices = cur_pos + moves_3d
#     # get edge locations
#     min_val = 0
#     max_val = grid.size  # not + 1 as in setup as we can walk on 0 or 100
#     # one will be at max and the other at min, that defines the type
#     for i in range(len(choices)):
#         temp = choices[i]
#         # 3d here, X reflective Y,Z periodic
#         if (temp[0] < min_val) or (temp[0] > max_val):  # X reflective, delete choice
#             del choices[i]
#             continue
#         #
#         if (temp[1] < min_val) and (temp[2] < min_val):
#             temp[1] = max_val
#             temp[2] = max_val
#         elif (temp[1] > max_val) and (temp[2] > max_val):
#             temp[1] = min_val
#             temp[2] = min_val
#         elif (temp[1] > max_val) and (temp[2] < min_val):
#             temp[1] = min_val
#             temp[2] = max_val
#         elif (temp[1] < min_val) and (temp[2] > max_val):
#             temp[1] = max_val
#             temp[2] = min_val
#         choices[i] = temp
#     return choices
""""""

"""Old kapitza cntvol 3D function where I tried to have detailed balance and P_m-cn ind. probabilities.
That doesn't actually make sense!"""

# def kapitza_cntvol(grid, moves_3d, kapitza, cur_pos, cur_index, prob_m_cn, inside_cnt):
#     random_num = np.random.random()  # [0.0, 1.0)
#     # check if the walker is inside or outside of a CNT
#     if inside_cnt:
#         # probs = [2.0 / 6.0, 4.0 / 6.0]  # detailed balance
#         kap_stay_enter = (random_num > prob_m_cn)
#         kap_leave_notenter = (random_num < prob_m_cn)
#     else:
#         #probs = [1.0 / 6.0, 5.0 / 6.0]  # detailed balance
#         kap_stay_enter = (random_num < prob_m_cn)
#         kap_leave_notenter = (random_num > prob_m_cn)
#     # d_b = ['stay_enter', 'leave_notenter']  # two possibilities
#     # d_b_choice = np.random.choice(d_b, p=probs)
#     if kap_stay_enter:
#         # (Detailed balance stay) AND (Kapitza stay)
#         # OR
#         # (Detailed balance leave) AND (Kapitza stay)
#         # move to random volume/endpoint within same CNT, remove current spot from choices
#         new_choices = []
#         for x in grid.tube_squares[cur_index - 1]:
#             if x not in [cur_pos]:
#                 new_choices.append(x)
#         num_new_choices = len(new_choices)
#         final_pos = np.asarray(new_choices[np.random.randint(0, num_new_choices)])
#     elif kap_leave_notenter:
#         # (Detailed balance stay) AND (Kapitza leave)
#         # OR
#         # (Detailed balance leave) AND (Kapitza leave)
#         # walk away, checking that current CNT volume is not a possibility
#         possible_locs, num_possible_locs = generate_novol_choices_3d(grid, moves_3d, cur_pos, cur_index, kapitza,
#                                                                      return_pos=True)
#         final_pos = np.asarray(possible_locs[np.random.randint(0, num_possible_locs)])
#     else:
#         kill()
#     #
#     inside_cnt = not inside_cnt
#     #
#     return final_pos, inside_cnt
""""""

"""Old apply moves code. This shows why the rules didnt make sense before 4/1/17. """
# def apply_moves_2d(walker, kapitza, grid, prob_m_cn, object):
#     '''Maybe given walker object or coordinates directly, check object'''
#     moves = [[0, 1], [1, 0], [0, -1], [-1, 0]]
#     # check where we are
#     if object:
#         cur_pos = np.asarray(walker.pos[-1])
#     else:
#         cur_pos = np.asarray(walker[-1])
#     # DEBUG
#     # Having tube radius doesn't matter if kapitza is off, apart from excluded volume
#     if kapitza:
#         cur_type = grid.tube_check_bd_vol[cur_pos[0], cur_pos[1]]  # type of square we're on
#         random_num = np.random.random()  # [0,1)
#         if cur_type == 1:  # endpoint
#             # step any direction
#             d_pos = np.asarray(moves[np.random.randint(0, 4)])
#             final_pos = cur_pos + d_pos
#             # walker.add_pos(final_pos)
#         elif cur_type == 0:  # matrix cell
#             # generate candidate position
#             d_pos = np.asarray(moves[np.random.randint(0, 4)])
#             candidate_pos = cur_pos + d_pos
#             candidate_type = grid.tube_check_bd_vol[candidate_pos[0], candidate_pos[1]]
#             if candidate_type == -1:  # inside tube
#                 if random_num > prob_m_cn:
#                     final_pos = cur_pos
#                     # walker.add_pos(final_pos)
#                 else:
#                     final_pos = candidate_pos
#                     # walker.add_pos(final_pos)
#             else:  # a normal step
#                 final_pos = candidate_pos
#                 # walker.add_pos(final_pos)
#         elif cur_type == -1:  # CNT cell
#             # find index of current tube walker is in
#             cur_index = grid.tube_check_index[cur_pos[0], cur_pos[1]]
#             cubes_in_tube = len(grid.tube_squares[cur_index - 1])
#             # -1 ABOVE AND BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
#             candidate_pos = grid.tube_squares[cur_index - 1][np.random.randint(0, cubes_in_tube)]
#             # move to another random point in the tube
#             d_pos = np.asarray(moves[np.random.randint(0, 4)])
#             candidate_pos = np.asarray(candidate_pos) + d_pos
#             # see where the candidate pos is
#             candidate_loc = grid.tube_check_bd_vol[candidate_pos[0], candidate_pos[1]]
#             # if candidate is in tube or on endpoint or random < kapitza move, else stay
#             if (candidate_loc == -1) or (candidate_loc == 1) or (random_num < prob_m_cn):
#                 final_pos = candidate_pos
#                 # walker.add_pos(final_pos)
#             else:
#                 final_pos = candidate_pos
#                 # walker.add_pos(final_pos)
#     else:  # standard hopping through tubes
#         cur_type = grid.tube_check_bd[cur_pos[0], cur_pos[1]]  # type of square we're on
#         if cur_type == 0:  # matrix cell
#             d_pos = np.asarray(moves[np.random.randint(0, 4)])
#             final_pos = cur_pos + d_pos
#             # walker.add_pos(final_pos)
#         elif cur_type == 1:  # endpoint
#             jump_moves = [[0, 1], [1, 0], [0, -1], [-1, 0], [0, 1], [1, 0], [0, -1], [-1, 0]]
#             # 8 possible jump positions for now if at tube end, 4 at left (first 4) and 4 at right (second 4)
#             choice = np.random.randint(0, 8)
#             d_pos = np.asarray(jump_moves[choice])
#             # coord on left tube end jumps to right end
#             tube_check_val_l = grid.tube_check_l[cur_pos[0], cur_pos[1]]
#             if tube_check_val_l > 0:
#                 if choice <= 3:  # stay at left end
#                     final_pos = cur_pos + d_pos
#                     # walker.add_pos(final_pos)
#                 else:  # jump across tube to right end
#                     # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
#                     final_pos = np.asarray(grid.tube_coords_r[tube_check_val_l - 1]) + np.asarray(d_pos)
#                     # walker.add_pos(final_pos)
#             # coord on right tube end jumps to left end
#             tube_check_val_r = grid.tube_check_r[cur_pos[0], cur_pos[1]]
#             if tube_check_val_r > 0:
#                 if choice <= 3:  # stay at right end
#                     final_pos = cur_pos + d_pos
#                     # walker.add_pos(final_pos)
#                 else:  # jump across tube to left end
#                     # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
#                     final_pos = np.asarray(grid.tube_coords_l[tube_check_val_r - 1]) + np.asarray(d_pos)
#                     # walker.add_pos(final_pos)
#     if cur_type == 10:  # reflective boundary
#         d_pos = np.asarray(moves[np.random.randint(0, 4)])
#         candidate_pos = cur_pos + d_pos
#         if (candidate_pos[0] > grid.size) or (candidate_pos[0] < 0):
#             final_pos = cur_pos
#             # walker.add_pos(final_pos)
#         else:
#             final_pos = candidate_pos
#             # walker.add_pos(final_pos)
#     elif cur_type == 20:  # periodic boundary
#         d_pos = np.asarray(moves[np.random.randint(0, 4)])
#         candidate_pos = cur_pos + d_pos
#         if candidate_pos[1] > grid.size:
#             final_pos = [candidate_pos[0], 1]
#             # walker.add_pos(final_pos)
#         elif candidate_pos[1] < 0:
#             final_pos = [candidate_pos[0], grid.size - 1]
#             # walker.add_pos(final_pos)
#         else:
#             final_pos = candidate_pos
#             # walker.add_pos(final_pos)
#     elif cur_type == 30:  # corner
#         d_pos = np.asarray(moves[np.random.randint(0, 4)])
#         candidate_pos = cur_pos + d_pos
#         if (candidate_pos[0] > grid.size) or (candidate_pos[0] < 0) or (candidate_pos[1] > grid.size) or (
#             candidate_pos[1] < 0):
#             final_pos = cur_pos
#             # walker.add_pos(final_pos)
#         else:
#             final_pos = candidate_pos
#     if object:
#         walker.add_pos(final_pos)
#     else:
#         walker.append(list(final_pos))
#     return walker
""""""
