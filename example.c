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
#define MAX_NUM_CELLS_PER_NET 500

typedef enum PARTITION {
    UNKNOWN,
    LEFT,
    RIGHT,
} PARTITION;

static char *PARTITION_STR[3] = { "UNKNOWN", "LEFT", "RIGHT" };

typedef struct LOGIC_CELL {
    int id;
    bool locked;
    int gain;
    PARTITION partition;

    struct LOGIC_CELL *next;
    struct LOGIC_CELL *prev;
} LOGIC_CELL;

LOGIC_CELL *logic_cells;

typedef struct NET {
    int num_cells;
    int net_id;

    LOGIC_CELL *cells[MAX_NUM_CELLS_PER_NET];

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
    INIT1,
    RECURSE,
    EXIT,
} STATE;

STATE state = IDLE;
bool done = false;
int num_iterations = -1;
int exit_counter = 0;

int num_cells = -1;
int balance = -1;
int num_cnx = -1;
int num_cols_per_partition = 0;
int num_rows_per_partition = 0;
int final_cost = 0;

void delay();
void button_press(float x, float y);
void proceed_state_button_func(void (*drawscreen_ptr) (void));
void debug_button_func(void (*drawscreen_ptr) (void));
void drawscreen();
int parse_file(char *file);
void print_net(NET *net);
void add_to_list(NET **head, NET *n);
void add_to_list(LOGIC_CELL **head, LOGIC_CELL *n);
LOGIC_CELL *make_logic_cell(int id, PARTITION p);
LOGIC_CELL *get_logic_cell(LOGIC_CELL *head, int id);
LOGIC_CELL *copy_logic_cells(LOGIC_CELL *l);
void destroy_logic_cell(LOGIC_CELL *l);

int calculate_cost(LOGIC_CELL *l);
int calculate_cost();

void init_grid();
double random(double from, double to);
float std_dev(int *vals, int size);
bool take_move(int delta_cost);

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

void reset_grid(PARTITION partition) {
    CELL **grid = (partition == LEFT) ? grid_left : grid_right;
    for (int col = 0; col < num_cols_per_partition; col++) {
        for (int row = 0; row < num_rows_per_partition; row++) {
            grid[col][row].cell_id = -1; 
        }
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

bool all_nodes_locked() {
    LOGIC_CELL *cur = logic_cells;
    while (cur != NULL) {
        if (!cur->locked) {
            return false;
        }
        cur = cur->next;
    }
    return true;
}

int get_gain(LOGIC_CELL *l) {
    // Gain of a node = # edges that cross partition - # edges that do not
    NET *cur = all_nets;

    int edges_cross = 0;
    int edges_dont_cross = 0;

    // Go through all nets
    while (cur != NULL) {
        PARTITION src_partition = UNKNOWN;

        // For each net, go through its connections
        for (int i = 0; i < cur->num_cells; i++) {

            // See if this net contains this cell
            if (l->id == cur->cells[i]->id) {
                if (i == 0) { // source cell
                    src_partition = l->partition;
                    continue;
                } else {
                    // Not a source...so find the source
                    LOGIC_CELL *c = get_logic_cell(logic_cells, cur->cells[0]->id);
                    src_partition = c->partition;

                    if (src_partition == l->partition) {
                        edges_dont_cross++;
                    } else {
                        edges_cross++;
                    }
                    break;
                }
            }

            // This is for if the cell is a source
            if (src_partition != UNKNOWN) {
                LOGIC_CELL *c = get_logic_cell(logic_cells, cur->cells[i]->id);
                if (c->partition == src_partition) {
                    edges_dont_cross++;
                } else {
                    edges_cross++;
                }
            }
        }

        cur = cur->next;
    }

    printf("Cell ID: %d has gain %d\n", l->id, edges_cross - edges_dont_cross);

    return (edges_cross - edges_dont_cross);
}

void init_partition() {
    int array[num_cells];
    for (int i = 0; i < num_cells; i++) {
        array[i] = i;
    }

    randomize_array(array, num_cells);

    for (int i = 0; i < num_cells; i++) {
        printf("array[%d]: %d\n", i, array[i]);
    }

    for (int i = 0; i < num_cells; i++) {
        LOGIC_CELL *cur = get_logic_cell(logic_cells, i);
        cur->partition = (array[i] % 2 == 0) ? LEFT : RIGHT;
    }
}

void init_partition2() {
    int best_cost = INT_MAX;

    LOGIC_CELL *tmp = copy_logic_cells(logic_cells);

    while (!all_nodes_locked()) {
        // Calculate all gains
        LOGIC_CELL *cur = tmp;
        int id_with_highest_gain = -1;
        int highest_gain = -INT_MAX;
        int num_right = 0;
        int num_left = 0;
        while (cur != NULL) {
            int gain = get_gain(cur);
            cur->gain = gain;
            if (cur->partition == LEFT) num_left++;
            else if (cur->partition == RIGHT) num_right++;
            cur = cur->next;
        }

        printf("Num left: %d Num right: %d\n", num_left, num_right);

        // Now go through and see which ones we can move
        cur = tmp;
        while (cur != NULL) {
            bool is_candidate = false;
            if (!cur->locked) {
                // Cell isn't locked...now see if we're still balanced
                if (num_left > num_right) {
                    // Then we can only move left side cells
                    if (cur->partition == LEFT) {
                        is_candidate = true;
                    }
                } else if (num_left < num_right) {
                    // Then we can only move right side cells
                    if (cur->partition == RIGHT) {
                        is_candidate = true;
                    }
                } else {
                    // Then who cares...
                    is_candidate = true;
                }
            }

            if (is_candidate) {
                printf("Cell %d is a candidate with gain %d vs highest gain: %d\n", cur->id, cur->gain, highest_gain);
                // If the candidate has a higher gain than what we've seen so far
                if (cur->gain > highest_gain) {
                    highest_gain = cur->gain;
                    id_with_highest_gain = cur->id;
                }
            }
            cur = cur->next;
        }

        // Choose node from candidate cells with highest gain
        if (id_with_highest_gain != -1) {
            LOGIC_CELL *c = get_logic_cell(tmp, id_with_highest_gain);
            if (c->partition == LEFT) {
                c->partition = RIGHT;
            } else {
                c->partition = LEFT;
            }

            int cost = calculate_cost(tmp);
            if (cost < best_cost) {
                best_cost = cost;
                // Record the solution
                LOGIC_CELL *l_cur = logic_cells;
                while (l_cur != NULL) {
                    LOGIC_CELL *c = get_logic_cell(tmp, l_cur->id);
                    l_cur->partition = c->partition;
                    l_cur = l_cur->next;
                }
            }

            c->locked = true;
        } else {
            printf("Can't find a candidate cell with the highest gain\n");
            break;
        }
    }
    destroy_logic_cell(tmp);
}

void reset_partition() {
    LOGIC_CELL *cur = logic_cells;
    while (cur != NULL) {
        cur->partition = UNKNOWN;
        cur = cur->next;
    }
}

int calculate_cost(LOGIC_CELL *l) {
    NET *cur = all_nets;
    int total_cost = 0;

    while (cur != NULL) {
        PARTITION cur_partition = UNKNOWN;

        for (int i = 0; i < cur->num_cells; i++) {
            int cur_id = cur->cells[i]->id;
            LOGIC_CELL *c = get_logic_cell(l, cur_id);
            if (c != NULL) {
                if (cur_partition == UNKNOWN) {
                    cur_partition = c->partition;
                } else {
                    // see if the cell has a valid partition that is different
                    if (c->partition != UNKNOWN) {
                        if (c->partition != cur_partition) {
                            total_cost++;
                            break;
                        }
                    }
                }
            }
        }

        cur = cur->next;
    }

#ifdef DEBUG
    printf("Total cost: %d\n", total_cost);
#endif

    return total_cost;
}

int calculate_cost() {
    /**
     * Cost is determined by the number of nets that cross the partition. If a
     * net connects cells that falls in both partitions, then that net contributes
     * a cost of 1 to the total crossing count. That is, the maximum cost is the
     * number of nets in this circuit.
     */
    NET *cur = all_nets;
    int total_cost = 0;

    while (cur != NULL) {
        PARTITION cur_partition = UNKNOWN;
        // Go through all the cells of this net
        for (int i = 0; i < cur->num_cells; i++) {
            LOGIC_CELL *c = cur->cells[i];
            // First time
            if (cur_partition == UNKNOWN) {
                cur_partition = c->partition;
            } else {
                // see if the cell has a valid partition that is
                // different
                if (c->partition != UNKNOWN) {
                    if (c->partition != cur_partition) {
                        total_cost++;
                        break;
                    }
                }
            }
        }
        cur = cur->next;
    }

#ifdef DEBUG
    printf("Total cost: %d\n", total_cost);
#endif
    return total_cost;
}

int size_of(LOGIC_CELL *head) {
    LOGIC_CELL *cur = head;
    int size = 0;
    while (cur != NULL) {
        size++;
        cur = cur->next;
    }
    return size;
}

bool is_valid_current_assignment(LOGIC_CELL *head) {
    bool is_valid = true;
    int num_left = 0;
    int num_right = 0;
    LOGIC_CELL *cur = head;
    while (cur != NULL) {
        if (cur->partition == LEFT) {
            num_left++;
        } else if (cur->partition == RIGHT) {
            num_right++;
        }

        if (num_left > balance || num_right > balance) {
            is_valid = false;
            break;
        }
        cur = cur->next;
    }
    return is_valid;
}

bool is_valid_solution(LOGIC_CELL *head) {
    bool is_valid = true;
    LOGIC_CELL *cur = head;
    if (size_of(head) == num_cells) {
        int num_left = 0;
        int num_right = 0;
        while (cur != NULL) {
            if (cur->partition == LEFT) {
                num_left++;
            } else if (cur->partition == RIGHT) {
                num_right++;
            } else {
                printf("WARNING: UNKNOWN partition!\n");
                return false;
            }

            if (num_left > balance || num_right > balance) {
                return false;
            }
            cur = cur->next;
        }
    } else {
        printf("Num cells: %d Size of cells: %d\n", num_cells, size_of(head));
        is_valid = false;
    }

    return is_valid;
}

void do_recursion(LOGIC_CELL *cur_ass, LOGIC_CELL *next_node) {

#ifdef DEBUG
    printf("Current assignment size: %d, Next node: %d Final_cost: %d\n", size_of(cur_ass), (next_node == NULL) ? -1 : next_node->id, final_cost);
#endif

    if (!is_valid_current_assignment(cur_ass)) {
        destroy_logic_cell(cur_ass);
        destroy_logic_cell(next_node);
#ifdef DEBUG
        printf("Pruning...not valid current assignment %d\n", final_cost);
#endif
        return;
    }
    // if there is no next node to assign
    if (next_node == NULL) {
        int cur_cost = calculate_cost(cur_ass);

        // if this is the best soln so far, record it
        if (cur_cost < final_cost) {
            final_cost = cur_cost;

            printf("Recorded better solution: %d\n", final_cost);

            // Record the partition solution
            LOGIC_CELL *cur = logic_cells;
            while (cur != NULL) {
                LOGIC_CELL *c = get_logic_cell(cur_ass, cur->id);
                cur->partition = c->partition;
                cur = cur->next;
            }
        }

        destroy_logic_cell(cur_ass);
        destroy_logic_cell(next_node);
    } else {
        // calculate label(x)
        int cur_cost = calculate_cost(cur_ass);

        // if (x < best solution so far) 
        if (cur_cost < final_cost) {
            int nn_id = next_node->id;
            LOGIC_CELL *nn_left = make_logic_cell(nn_id, LEFT);
            LOGIC_CELL *cur_left = copy_logic_cells(cur_ass);
            bool last_node = (nn_id + 1 >= num_cells);

            destroy_logic_cell(cur_ass);
            destroy_logic_cell(next_node);
            
            LOGIC_CELL *cur_right = copy_logic_cells(cur_left);
            LOGIC_CELL *nn_right = make_logic_cell(nn_id, RIGHT);

            add_to_list(&cur_left, nn_left);

            LOGIC_CELL *tmp_nn_left = NULL;
            if (!last_node) {
                tmp_nn_left = make_logic_cell(nn_id + 1, UNKNOWN);
            }
            do_recursion(cur_left, tmp_nn_left);

            add_to_list(&cur_right, nn_right);
            LOGIC_CELL *tmp_nn_right = NULL;
            if (!last_node) {
                tmp_nn_right = make_logic_cell(nn_id + 1, UNKNOWN);
            }
            do_recursion(cur_right, tmp_nn_right);
        } else {
            destroy_logic_cell(cur_ass);
            destroy_logic_cell(next_node);
        }
    }
}

void run_algo() {
    LOGIC_CELL *l = make_logic_cell(0, UNKNOWN);
    do_recursion(NULL, l);
}

void run_partition() {
    switch (state) {
        case IDLE: {
            printf("State is IDLE\n");
            init_partition();
            LOGIC_CELL *cur = logic_cells;
            while (cur != NULL) {
                assign_grid(cur->id, cur->partition);
                cur = cur->next;
            }
            final_cost = calculate_cost();
            state = INIT1;
        } break;
        case INIT1: {
            printf("State is INIT1\n");
            reset_grid(LEFT);
            reset_grid(RIGHT);

            init_partition2();
            LOGIC_CELL *cur = logic_cells;
            while (cur != NULL) {
                assign_grid(cur->id, cur->partition);
                cur = cur->next;
            }
            final_cost = calculate_cost();
            state = RECURSE;
        } break;
        case RECURSE: {
            printf("State is RECURSE\n");
            run_algo();
            reset_grid(LEFT);
            reset_grid(RIGHT);
            LOGIC_CELL *cur = logic_cells;
            while (cur != NULL) {
                assign_grid(cur->id, cur->partition);
                cur = cur->next;
            }
            state = EXIT;
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

void print_net(NET *net) {
    int i;
    if (net == NULL) {
        printf("Net is NULL\n");
    } else {
        printf("Net: num logic blocks: %d\n", net->num_cells);
        for (i = 0; i < net->num_cells; i++) {
            printf("     Cell[%d]: %d partition (%s)\n", i, net->cells[i]->id, PARTITION_STR[net->cells[i]->partition]);
        }
    }
}


void debug_button_func(void (*drawscreen_ptr) (void)) {
    printf("Logic cells:\n");
    LOGIC_CELL *tmp = logic_cells;
    while (tmp != NULL) {
        printf("  id: %d - partition (%s)\n", tmp->id, PARTITION_STR[tmp->partition]);
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
                printf("Parse_file[%d] (%d): %s", line_num, (int)strlen(line), line);
#endif

                // First line contains # cells, # cnx btw cells, # rows, # cols
                if (line_num == 0) {
                    const char delim[2] = " ";
                    char *token;

                    token = strtok(line, delim);
                    num_cells = atoi(token);
                    balance = num_cells / 2;
                    if (num_cells % 2) {
                        balance++;
                    }
                    token = strtok(NULL, delim);
                    num_cnx = atoi(token);

                    for (int i = 0; i < num_cells; i++) {
                        LOGIC_CELL *c = make_logic_cell(i, UNKNOWN);
                        add_to_list(&logic_cells, c);
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
                    net->num_cells = atoi(token);
                    net->net_id = line_num - 1;

                    // Remaining numbers are the block numbers connected to this net.
                    for (i = 0; i < net->num_cells; i++) {
                        token = strtok(NULL, delim);
                        int cell_id = atoi(token);
                        LOGIC_CELL *c = get_logic_cell(logic_cells, cell_id);
                        net->cells[i] = c;
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
    for (i = 0; i < 100; i++)
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

void add_to_list(LOGIC_CELL **head, LOGIC_CELL *n) {
    LOGIC_CELL *cur = *head;

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

LOGIC_CELL *make_logic_cell(int id, PARTITION p) {
    LOGIC_CELL *l = (LOGIC_CELL *)malloc(sizeof(LOGIC_CELL));
    l->id = id;
    l->partition = p;
    l->gain = 0;
    l->locked = false;
    l->next = NULL;
    l->prev = NULL;

    return l;
}

LOGIC_CELL *get_logic_cell(LOGIC_CELL *head, int id) {
    LOGIC_CELL *cur = head;
    LOGIC_CELL *ret = NULL;
    while (cur != NULL) {
        if (cur->id == id) {
            ret = cur;
            break;
        }
        cur = cur->next;
    }
    return ret;
}

LOGIC_CELL *copy_logic_cells(LOGIC_CELL *l) {
    LOGIC_CELL *cur = l;
    LOGIC_CELL *copy = NULL;
    while (cur != NULL) {
        LOGIC_CELL *tmp = make_logic_cell(cur->id, cur->partition);
        add_to_list(&copy, tmp);
        cur = cur->next;
    }

    return copy;
}

void destroy_logic_cell(LOGIC_CELL *l) {
    LOGIC_CELL *cur = l;
    while (cur != NULL) {
        LOGIC_CELL *next = cur->next;
        free(cur);
        cur = next;
    }
}

void button_press(float x, float y) {
    printf("User clicked a button at coordinates (%f, %f)\n", x, y);
}
