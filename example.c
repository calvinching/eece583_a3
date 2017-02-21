#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include "graphics.h"
#include <sys/time.h>

#define DEBUG
#define SUCCESS 0
#define ERROR -1
#define MAX_NUM_SINKS 500
#define NUM_INIT_ITERATIONS 50 // Number of times to run the random move to find the initial temperature
#define BETA_TEMP_FACTOR 0.9
#define EXIT_COUNT 500

#define IS_MIN_AND_SET(x, y) if (x < y) y = x;
#define IS_MAX_AND_SET(x, y) if (x > y) y = x;

typedef enum PARTITION {
    UNKNOWN,
    LEFT,
    RIGHT,
} PARTITION;

typedef struct LOGIC_CELL {
    int col;
    int row;

    int id;
    PARTITION partition;
} LOGIC_CELL;

LOGIC_CELL *logic_cells;

typedef struct NET {
    int num_logic_blocks;
    int net_id;

    LOGIC_CELL *src_cell;
    LOGIC_CELL *sink_cell[MAX_NUM_SINKS];

    struct NET *next;
    struct NET *prev;
} NET;

NET *all_nets = NULL;

typedef struct CELL {
    float x1;       // x-coordinate of the cell's top left corner
    float y1;       // y-coordinate of the cell's top left corner
    float x2;       // x-coordinate of the cell's bottom right corner
    float y2;       // y-coordinate of the cell's bottom right corner
    char text[5];   // text of the cell
    float text_x;   // x-coordinate of the text
    float text_y;   // y-coordinate of the text

    int cell_id;
} CELL;

CELL **grid_left;
CELL **grid_right;

typedef enum STATE {
    IDLE,
    INIT,
    ITERATE,
    PENDING_EXIT,
    EXIT,
} STATE;

STATE state = IDLE;
bool done = false;
int num_iterations = -1;
int exit_counter = 0;

int num_cells = -1;
int num_cnx = -1;
int num_cols_per_partition = 0;
int num_rows_per_partition = 0;

float temp = 0.;

void delay();
void button_press(float x, float y);
void proceed_button_func(void (*drawscreen_ptr) (void));
void proceed_state_button_func(void (*drawscreen_ptr) (void));
void debug_button_func(void (*drawscreen_ptr) (void));
void drawscreen();
int parse_file(char *file);
void print_net(NET *net);
void add_to_list(NET **head, NET *n);
void init_grid();
double random(double from, double to);
float std_dev(int *vals, int size);
bool take_move(int delta_cost);
void find_random_cells(int *c1, int *c2);


void init_grid(bool left) {
    t_report report;
    report_structure(&report);
#ifdef DEBUG
    printf("width: %d, height: %d\n", report.top_width, report.top_height);
#endif

    float height = report.top_height;
    float width = report.top_width;
    float cell_height = height / num_rows_per_partition;
    float cell_width = width / (num_cols_per_partition * 2);

    int offset = (left) ? 0 : width / 2 + 10;
    CELL **grid;

#ifdef DEBUG
    printf("offset: %d\n", offset);
#endif

    // Allocate memory for the grid
    grid = (CELL **)malloc(num_cols_per_partition * sizeof(CELL *));
    for (int col = 0; col < num_cols_per_partition; col++) {
        grid[col] = (CELL *)malloc(num_rows_per_partition * sizeof(CELL));
    }

    for (int col = 0; col < num_cols_per_partition; col++) {
        for (int row = 0; row < num_rows_per_partition; row++) {
            grid[col][row].x1 = cell_width * col + offset;
            grid[col][row].y1 = cell_height * row;
            grid[col][row].x2 = cell_width + cell_width * col + offset;
            grid[col][row].y2 = cell_height + cell_height * row;
            grid[col][row].text_x = grid[col][row].x2 - cell_width / 2.;
            grid[col][row].text_y = grid[col][row].y2 - cell_height / 2.;

            grid[col][row].cell_id = -1;
#ifdef DEBUG
            printf("grid_%s[%d][%d] = (%f, %f) (%f, %f) (%f, %f)\n", ((left) ? "left" : "right"), col, row, grid[col][row].x1, grid[col][row].y1, grid[col][row].x2, grid[col][row].y2, grid[col][row].text_x, grid[col][row].text_y);
#endif
        }
    }

    if (left) {
        grid_left = grid;
    } else {
        grid_right = grid;
    }
}

void init_grid() {
    int num_cells_per_partition = num_cells / 2;
    num_cols_per_partition = 0;
    num_rows_per_partition = 0;

    if (num_cells % 2) {
        num_cells_per_partition++;
    }

    for (int i = num_cells_per_partition; i > num_cells_per_partition / 2 - 1; i--) {
        int mod = num_cells_per_partition % i;
        if (mod == 0) {
            num_rows_per_partition = i;
        }
    }

    num_cols_per_partition = num_cells_per_partition / num_rows_per_partition;

    printf("Each partition: %d %d\n", num_cols_per_partition, num_rows_per_partition);
    init_grid(true);
    init_grid(false);
}

int main(int argc, char **argv) {
    if (argc != 2) {
        printf("Need input file\n");
        exit(1);
    }
    char *file = argv[1];
    printf("Input file: %s\n", file);

    // initialize display
    init_graphics("Some Example Graphics");
    init_world(0.,0.,1000.,1000.);

    parse_file(file);
    init_grid();

    create_button("Window", "Go 1 Step", proceed_button_func);
    create_button("Window", "Go 1 State", proceed_state_button_func);
    create_button("Window", "Debug", debug_button_func);
    drawscreen();
    event_loop(button_press, drawscreen);
    return 0;
}


double random(double from, double to) {
    struct timeval t1;
    gettimeofday(&t1, NULL);
    srand(t1.tv_usec * t1.tv_sec);

    return (rand() / (double)(RAND_MAX)) * abs(from - to) + from;
}

void randomize_array(int array[], int size) {
    struct timeval t1;
    gettimeofday(&t1, NULL);
    srand(t1.tv_usec * t1.tv_sec);

    for (int i = size - 1; i > 0; i--) {
        int j = rand() % (i+1);
        int tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
}

void assign_grid(int cell_id, PARTITION partition) {
    CELL **grid = (partition == LEFT) ? grid_left : grid_right;
    for (int col = 0; col < num_cols_per_partition; col++) {
        for (int row = 0; row < num_rows_per_partition; row++) {
            if (grid[col][row].cell_id == -1) {
                grid[col][row].cell_id = cell_id;
                printf("Assigned grid_%s[%d][%d] for cell_id: %d\n", (partition == LEFT) ? "left" : "right", col, row, cell_id);
                return;
            }
        }
    }
}

void init_partition() {
    for (int i = 0; i < num_cells; i++) {
        logic_cells[i].partition = (i % 2 == 0) ? LEFT : RIGHT;
        assign_grid(i, logic_cells[i].partition);
    }
    state = INIT;
/*
    int total = num_cols * num_rows;
    int array[total];
    for (int i = 0; i < total; i++) {
        array[i] = i;
    }

    randomize_array(array, total);

    for (int i = 0; i < num_cells; i++) {
        printf("array[%d]: %d\n", i, array[i]);
        int col = array[i] / num_rows;
        int row = array[i] % num_rows;

        if (grid[col][row].is_block) {
            printf("ERROR!!! (%d, %d)\n", col, row);
        } else {
            grid[col][row].is_block = true;
            grid[col][row].block_num = i;

            logic_cells[i].col = col;
            logic_cells[i].row = row;
            printf("Found (%d, %d) for logic block: %d\n", col, row, i);
        }
    }

    state = INIT;
    */
}

void find_random_cells(int *c1, int *c2) {
    int array[num_cells];
    for (int i = 0; i < num_cells; i++) {
        array[i] = i;
    }

    randomize_array(array, num_cells);

    *c1 = array[0];
    *c2 = array[1];
#ifdef DEBUG
    printf("Found random logic cell: %d, %d\n", *c1, *c2);
#endif
}

void swap_cells(LOGIC_CELL *c1, LOGIC_CELL *c2) {
    int c1_col = c1->col;
    int c1_row = c1->row;
    int c2_col = c2->col;
    int c2_row = c2->row;

/*
    int c1_blk = grid[c1_col][c1_row].block_num;
    int c2_blk = grid[c2_col][c2_row].block_num;
    grid[c1_col][c1_row].block_num = c2_blk;
    grid[c2_col][c2_row].block_num = c1_blk;

    c1->col = c2_col;
    c1->row = c2_row;
    c2->col = c1_col;
    c2->row = c1_row;
*/
#ifdef DEBUG
    printf("Swapped LOGIC_CELL %d from (%d, %d) to (%d, %d)\n", c1->id, c1_col, c1_row, c1->col, c1->row);
    printf("Swapped LOGIC_CELL %d from (%d, %d) to (%d, %d)\n", c2->id, c2_col, c2_row, c2->col, c2->row);
#endif
}


int calculate_cost() {
    /**
     * Calcuate cost:
     *
     * So suppose you have n nets, each with some number of pins (at least 2, but sometimes more).
     * For net i, we would cycle through all the pins in this net. Each pin has a location, xj and yj.  Suppose
     * the minimum xj for all pins on net i is xi,min and the maximum xj for all pins on net i is xi,max.
     * Further, suppose the minimum yj for all pins on net i is yi,min and the maximum yj for all pins on net i is yi,max.
     * Then, the size of net i can be approximated by
     *      si = (yi,max - yi,min ) + ( xi,max - xi,min )
     *
     * In other words, for each net, find the minimum rectangular area that encompasses all pins in that net.
     * Then si is the sum of the x span plus the y span of this rectangular area.  It is also 1/2 the perimeter of
     * the rectangular area (which is why it is called .half perimeter.).
     *
     * Then, the total cost of a placement is
     *      the sum of si for all i
     *          (in other words, add up the costs of each net to get an overall cost).
     */
    NET *cur = all_nets;
    int total_cost = 0;

    while (cur != NULL) {
        int cur_cost = 0;
        int col_min = INT_MAX;
        int col_max = 0;
        int row_min = INT_MAX;
        int row_max = 0;

        // Start with the source cell of this net
        if (cur->src_cell != NULL) {
            LOGIC_CELL *c = cur->src_cell;
            int col = c->col;
            int row = c->row;

            // Update min/max values for columns and rows
            IS_MIN_AND_SET(col, col_min);
            IS_MAX_AND_SET(col, col_max);
            IS_MIN_AND_SET(row, row_min);
            IS_MAX_AND_SET(row, row_max);
        } else {
            printf("WARN: invalid source logic cell\n");
        }

        // Go through all the sink cells of this net
        for (int i = 0; i < cur->num_logic_blocks - 1; i++) {
            LOGIC_CELL *c = cur->sink_cell[i];
            int col = c->col;
            int row = c->row;

            // Update min/max values for columns and rows
            IS_MIN_AND_SET(col, col_min);
            IS_MAX_AND_SET(col, col_max);
            IS_MIN_AND_SET(row, row_min);
            IS_MAX_AND_SET(row, row_max);
        }

        cur_cost = (col_max - col_min) + (row_max - row_min);
        total_cost += cur_cost;
#ifdef DEBUG
//        printf("Cost of net %d: (%d - %d) + (%d - %d) = %d\n", cur->net_id, col_max, col_min, row_max, row_min, cur_cost);
//        printf("Total cost: %d\n", total_cost);
#endif
        cur = cur->next;
    }
    return total_cost;
}

float std_dev(int *vals, int size) {
    float mean, stddev;
    float sum = 0;
    for (int i = 0; i < size; i++) {
        sum += (float)vals[i];
    }
    mean = sum / (float)size;
#ifdef DEBUG
    printf("Mean: %f\n", mean);
#endif
    sum = 0;
    for (int i = 0; i < size; i++) {
       sum += pow((float)vals[i] - mean, 2);
    }
    mean = sum / (float)size;
    stddev = sqrt(mean);
#ifdef DEBUG
    printf("Std dev: %f\n", stddev);
#endif
    return stddev;
}

void run_init_temp() {
    int c1 = -1;
    int c2 = -1;

    int all_costs[NUM_INIT_ITERATIONS];

    for (int i = 0; i < NUM_INIT_ITERATIONS; i++) {
        printf("Iteration: %d\n", i);

        // Find 2 random cells to swap
        find_random_cells(&c1, &c2);
        LOGIC_CELL *cell1 = &logic_cells[c1];
        LOGIC_CELL *cell2 = &logic_cells[c2];

        // Swap the two cells
        swap_cells(cell1, cell2);

        // Calculate the cost
        all_costs[i] = calculate_cost();
    }

    // Initial temp = 20 * std_dev(costs over the NUM_INIT_ITERATIONS moves
    temp = 20 * std_dev(all_costs, NUM_INIT_ITERATIONS);

    printf("Initial temperature: %f\n", temp);

    state = ITERATE;
}

void update_temp() {
    printf("Temperature: %f\n", temp);
    temp = temp * BETA_TEMP_FACTOR;
    printf("New temperature: %f\n", temp);
}

void run_partition() {
    switch (state) {
        case IDLE: {
            printf("State is IDLE\n");
            init_partition();
        } break;
        case INIT: {
            printf("State is INIT\n");
            run_init_temp();
            num_iterations = 10 * pow((double)num_cells, (double)(4.0/3.0));
            printf("Number of iterations: %d\n", num_iterations);
        } break;
        case ITERATE: {
            int cur_cost, prev_cost = -1;
            for (int i = 0; i < num_iterations; i++) {
                int c1, c2;
                find_random_cells(&c1, &c2);
                LOGIC_CELL *cell1 = &logic_cells[c1];
                LOGIC_CELL *cell2 = &logic_cells[c2];

                // Get previous cost
                if (prev_cost == -1) {
                    prev_cost = calculate_cost();
                }

                // Swap the two cells
                swap_cells(cell1, cell2);

                // Get the cost of the new placement
                cur_cost = calculate_cost();

                // Find the change in cost
                int delta_cost = cur_cost - prev_cost;
                printf("cur_cost: %d - prev_cost: %d = delta_cost: %d\n", cur_cost, prev_cost, delta_cost);
                if (take_move(delta_cost)) {
                    printf("Taking move in iteration %d\n", i);

                    // Update cost
                    prev_cost = cur_cost;

                    // Reset the exit counter if we are taking a move
                    exit_counter = 0;
                } else {
                    printf("***Not taking the move in iteration %d exit_counter: %d temp: %f\n", i, exit_counter, temp);
                    // Undo the swap
                    swap_cells(cell1, cell2);
                    exit_counter++;
                }

                // Check to see if we need to update state (i.e. meet exit criteria)
                if (exit_counter > EXIT_COUNT) {
                    printf("We are exiting due to exit_counter. Setting state to PENDING_EXIT\n");
                    printf("Cost: %d Temp: %f\n", prev_cost, temp);
                    state = PENDING_EXIT;
                    break;
                }
            }

            if (state != PENDING_EXIT) {
                printf("Finished one iteration loop for temp: %f. Cost: %d\n", temp, cur_cost);
                // Update temp
                update_temp();
            }
        } break;
        case PENDING_EXIT: {
            // Before we quit, let's make one more swap and see if we improve
            int c1, c2;
            find_random_cells(&c1, &c2);
            LOGIC_CELL *cell1 = &logic_cells[c1];
            LOGIC_CELL *cell2 = &logic_cells[c2];

            // Get previous cost
            int prev_cost = calculate_cost();

            // Swap the two cells
            swap_cells(cell1, cell2);

            // Get the cost of the new placement
            int cur_cost = calculate_cost();

            // Find the change in cost
            int delta_cost = cur_cost - prev_cost;

            printf("%d - %d = delta_cost %d\n", cur_cost, prev_cost, delta_cost);

            if (delta_cost < 0) {
                // We are still improving so we should keep iterating
                printf("Resetting state to ITERATE\n");
                state = ITERATE;
            } else {
                state = EXIT;
                printf("Final cost: %d Final Temp: %f\n", prev_cost, temp);
                debug_button_func(NULL);
            }
        } break;
        case EXIT: {
            printf("We are done!\n");
            done = true;
        } break;
        default: {
            printf("ERROR: unknown state!\n");
        } break;
    }
}

void proceed_button_func(void (*drawscreen_ptr) (void)) {
    run_partition();
    drawscreen();
}

void proceed_state_button_func(void (*drawscreen_ptr) (void)) {
    STATE cur_state = state;
    while (cur_state == state && !done) {
        run_partition();
        drawscreen();
        delay();
    }

    if (done)
        printf("Nothing else to do!\n");
}


bool take_move(int delta_cost) {
    // If cost is less than 0, then we take move for sure
    if (delta_cost < 0) {
        return true;
    }
    // If no cost change, just ignore it.
    if (delta_cost == 0) {
        return false;
    }

    // Get a random number between 0 and 1
    double rand = random(0, 1);

    // We take move if rand < e^(-delta_cost/temp)
    double e = exp((double)((delta_cost * -1.0)/(double)temp));
    printf("rand %f < exp: %f ?\n", rand, e);
    return (rand < e);
}

void print_net(NET *net) {
    int i;
    if (net == NULL) {
        printf("Net is NULL\n");
    } else {
        printf("Net: num logic blocks: %d\n", net->num_logic_blocks);
        printf("     Source: %d (%d, %d)\n", net->src_cell->id, net->src_cell->col, net->src_cell->row);
        for (i = 0; i < net->num_logic_blocks - 1; i++) {
            printf("     Sink[%d]: %d (%d, %d)\n", i, net->sink_cell[i]->id, net->sink_cell[i]->col, net->sink_cell[i]->row);
        }
    }
}


void debug_button_func(void (*drawscreen_ptr) (void)) {
    printf("Logic cells:\n");
    for (int i = 0; i < num_cells; i++) {
        printf("  [%d] - id: %d (%d, %d)\n", i, logic_cells[i].id, logic_cells[i].col, logic_cells[i].row);
    }

    printf("\nNets: \n");
    NET *cur = all_nets;
    while (cur != NULL) {
        print_net(cur);
        cur = cur->next;
    }
}

int parse_file(char *file) {
    FILE *fp;
    int ret = SUCCESS;

    if (file != NULL) {
        fp = fopen(file, "r");
        if (fp == NULL) {
            printf("Failed to open file: %s\n", file);
            ret = ERROR;
        } else {
            char *line;
            size_t len = 0;
            ssize_t read;

            int line_num = 0;
            while ((read = getline(&line, &len, fp)) != -1) {
                if (strlen(line) < 4) {
                    continue;
                }
#ifdef DEBUG
                printf("Parse_file[%d] (%d): %s", line_num, strlen(line), line);
#endif

                // First line contains # cells, # cnx btw cells, # rows, # cols
                if (line_num == 0) {
                    const char delim[2] = " ";
                    char *token;

                    token = strtok(line, delim);
                    num_cells = atoi(token);
                    token = strtok(NULL, delim);
                    num_cnx = atoi(token);

                    logic_cells = (LOGIC_CELL *)malloc(sizeof(LOGIC_CELL) * num_cells);
                    for (int i = 0; i < num_cells; i++) {
                        logic_cells[i].id = i;
                        logic_cells[i].col = -1;
                        logic_cells[i].row = -1;
                        logic_cells[i].partition = UNKNOWN;
                    }
#ifdef DEBUG
                    printf("Num cells: %d, Num cnx: %d\n", num_cells, num_cnx);
#endif
                } else {
                    // Remaining lines indicate nets. Each net can connect to 2 or more logic blocks.
                    const char delim[2] = " ";
                    char *token;
                    int i;
                    NET *net = (NET *)malloc(sizeof(NET));
                    net->next = NULL;
                    net->prev = NULL;

                    // First number is # logic blocks this net connects
                    token = strtok(line, delim);
                    net->num_logic_blocks = atoi(token);
                    net->net_id = line_num - 1;

                    // Remaining numbers are the block numbers connected to this net.
                    // Get the source block first
                    token = strtok(NULL, delim);
                    int src_id = atoi(token);
                    net->src_cell = &(logic_cells[src_id]);

                    // Rest are sink blocks
                    for (i = 0; i < net->num_logic_blocks - 1; i++) {
                        token = strtok(NULL, delim);
                        int snk_id = atoi(token);
                        net->sink_cell[i] = &(logic_cells[snk_id]);
                    }

#ifdef DEBUG
                    print_net(net);
#endif
                    add_to_list(&all_nets, net);
                }
                line_num++;
            }
        }
    } else {
        printf("Invalid file!\n");
    }

    return ret;
}

void draw_grid(bool left) {
    CELL **grid = (left) ? grid_left : grid_right;
    for (int col = 0; col < num_cols_per_partition; col++) {
        for (int row = 0; row < num_rows_per_partition; row++) {
            char text[10] = "";
            if (grid[col][row].cell_id != -1) {
                // Draw source and sinks
                setcolor(BLACK);
                drawrect(grid[col][row].x1, grid[col][row].y1, grid[col][row].x2, grid[col][row].y2);
                sprintf(text, "%d", grid[col][row].cell_id);
                drawtext(grid[col][row].text_x, grid[col][row].text_y, text, 150.);
            } else {
                setcolor(BLACK);
                drawrect(grid[col][row].x1, grid[col][row].y1, grid[col][row].x2, grid[col][row].y2);
#ifdef DEBUG
                sprintf(text, "%d,%d", col, row);
                drawtext(grid[col][row].text_x, grid[col][row].text_y, text, 150.);
#endif
            }
        }
    }
}

void drawscreen() {
    clearscreen();
    draw_grid(true);
    draw_grid(false);
}

void delay() {
    int i, j, k, sum;

    sum = 0;
    for (i = 0; i < 1000; i++)
        for (j = 0; j < i; j++)
            for (k = 0; k < 30; k++)
                sum = sum + i + j - k;
}

void add_to_list(NET **head, NET *n) {
    NET *cur = *head;

    if (*head == NULL) {
        *head = n;
        return;
    }

    // Go to the end of the list
    while (cur->next != NULL) {
        cur = cur->next;
    }

    cur->next = n;
    n->prev = cur;
    n->next = NULL;
}


void button_press(float x, float y) {
    printf("User clicked a button at coordinates (%f, %f)\n", x, y);
}
