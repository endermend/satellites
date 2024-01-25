#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <chrono>

#define earthRadius 6378135.0
#define Alpha 0.73303828583761842230795012276522

using namespace std;

//Клетки
typedef struct _Cell
{
    double x, y, z;
    double theta, phi;
    int sattelite_id;
    int max_satellite;
    int connection_time;
} CELL,*HCELL;

//Спутник
typedef struct _Satellite
{
    double x, y, z;
    double phi, theta;
    double heigth, delta;
} SATELLITE, *HSATELLITE;

//Сектора клеток
typedef struct _Sector
{
    HCELL *cells;
    int cnt;
    int gap;
} SECTOR,*HSECTOR;

HSECTOR grid = NULL;//Сетка клеток
HCELL cells = NULL;//Клетки
HSATELLITE satellites = NULL;//Спутники

int cellsTotalCount = 0;//Общее колличество клеток
int satellitesTotalCount = 0;//Колличество спутников
int timesTotalCount = 0;//Колличество отсчетов времени


//Добавление новой клетки в сектор
int push_to_Sector(HSECTOR sector, HCELL cell)
{
    if(!sector || !cell)
    {
        cout<<"Can't add cell" <<endl;
        return 0;
    }
    if(!sector->cells || sector->gap == 0)
    {
        sector->gap = 4;
        sector->cells = (HCELL *)malloc(sector->gap * sizeof(HCELL));
    }
    if(sector->cnt == sector->gap)
    {
        sector->gap <<= 1;
        sector->cells = (HCELL *)realloc(sector->cells,sector->gap*sizeof(HCELL));
    }
    sector->cells[sector->cnt++] = cell;
    return 1;
}

//Перевод системы координат клетки в сферическую и добавление клетки в сетку
int proc_Cell(HCELL cell, HSECTOR grid)
{
    if(!grid || !cell)
    {
        cout<<"Can't process cell" <<endl;
        return 0;
    }

    cell->sattelite_id = -1;
    cell->max_satellite = -1;
    cell->connection_time = 0;

    //Перевод координат в сферические
    cell->phi = atan(cell->y / cell->x);
    cell->theta = acos(cell->z / earthRadius);

    //Распределение угла по всему кругу
    if (cell->x > 0 && cell->y >= 0);
    else if (cell->x <= 0)
        cell->phi += M_PI;
    else if (cell->x >= 0 && cell->y < 0)
        cell->phi += 2 * M_PI;

    //Перевод радиан в градусы, для удобства представления
    cell->phi *= 180/M_PI;
    cell->theta *= 180/M_PI;

    //Добавление клетки в сетку
    push_to_Sector(grid+((int)(cell->phi))*180+(int)(cell->theta),cell);
    return 1;
}

//Перевод системы координат спутника в сферическую и расчет дальности обзора
int proc_Satellite(HSATELLITE satellite)
{
    if(!satellite)
    {
        cout << "Can't process satellite" << endl;
        return 0;
    }

    //Перевод координат в сферические
    satellite->phi = atan(satellite->y / satellite->x);
    satellite->heigth = sqrt(satellite->x*satellite->x + satellite->y*satellite->y + satellite->z*satellite->z);
    satellite->theta = acos(satellite->z / satellite->heigth);

    //Расчет дальности обзора delta
    double _alpha = M_PI / 2.0 - Alpha;
    double _cos2 = cos(_alpha) * cos(_alpha);
    double _x = (sqrt(earthRadius * earthRadius / _cos2 - satellite->heigth * satellite->heigth) - satellite->heigth * tan(_alpha)) * _cos2;
    double _z = tan(_alpha) * _x + satellite->heigth;
    satellite->delta = abs(atan(_x / _z));

    //Распределение угла по всему кругу
    if (satellite->x > 0 && satellite->y >= 0);
    else if (satellite->x <= 0)
        satellite->phi += M_PI;
    else if (satellite->x >= 0 && satellite->y < 0)
        satellite->phi += 2 * M_PI;

    //Перевод радиан в градусы, для удобства представления
    satellite->phi*=180.0/M_PI;
    satellite->theta*=180.0/M_PI;
    satellite->delta*=180.0/M_PI;

    return 1;
}

//Возвращает угол между спутником и клеткой относительно центра Земли
double calc_angle(HSATELLITE satellite, HCELL cell)
{
    double dot = satellite->x * cell->x + satellite->y * cell->y + satellite->z * cell->z;
    return abs(acos(dot / satellite->heigth / earthRadius))*180.0/M_PI;
}

//Расчитывает, сколько отсчётов времени точка еще потдерживается
int predict_support_time(HCELL cell, HSATELLITE satellite, HSATELLITE satellites_end, int satellitesTotalCount){
    HSATELLITE future_sat = satellite+satellitesTotalCount*128;
    for(; future_sat<=satellites_end && calc_angle(future_sat,cell) <= future_sat->delta; future_sat+=satellitesTotalCount*128);
    future_sat -= satellitesTotalCount*128;
    for(int delta = 64; delta > 0; delta>>=1)
    {
        if(future_sat+satellitesTotalCount*delta <= satellites_end && calc_angle(future_sat+satellitesTotalCount*delta,cell) <= (future_sat+satellitesTotalCount*delta)->delta)
        {
            future_sat+=satellitesTotalCount*delta;
        }
        if((future_sat-satellite)/satellitesTotalCount+delta-1 <= cell->connection_time)
        {
            break;
        }
    }
    return (future_sat-satellite)/satellitesTotalCount;
}

//Проверяет каждую клетку в секторе на покрытие ее спутником
int proc_sector(HSECTOR sector, HSATELLITE satellite, HSATELLITE satellites_end, int satellitesTotalCount)
{
    if(!sector || !sector->cells || !satellite)
    {
        cout<<"Can't process sector"<<endl;
        return 0;
    }
    for(HCELL *cell = sector->cells; cell < sector->cells+sector->cnt; cell++)
    {
        if((*cell)->sattelite_id > -1)
        {
            continue;
        }
        if(calc_angle(satellite,(*cell)) > satellite->delta)
        {
            continue;
        }
        //Если время которое спутник потдерживает точку больше, чем у другого спутника, то клетка выбирает этот спутник
        int observed_times = predict_support_time(*cell,satellite,satellites_end,satellitesTotalCount);
        if(observed_times>(*cell)->connection_time)
        {
            (*cell)->connection_time = observed_times;
            (*cell)->max_satellite = satellitesTotalCount-((satellites_end-satellite)%satellitesTotalCount);
        }
    }
}

int proc_time(HSATELLITE current_satellite, HSATELLITE satellites_end, int satellitesTotalCount, HSECTOR grid)
{
    if(!current_satellite || !satellites_end || !grid || current_satellite > satellites_end)
    {
        cout<<"Can't process satellite calculation"<<endl;
        return 0;
    }
    int lower_phi = (int)(current_satellite->phi-current_satellite->delta);
    int upper_phi = (int)(current_satellite->phi+current_satellite->delta)+1;

    int lower_theta = (int)(current_satellite->theta-current_satellite->delta);
    int upper_theta = (int)(current_satellite->theta+current_satellite->delta)+1;

    //случай если спутник покрывает северный полюс
    if(upper_theta >180){
        for(int cell_phi = 0; cell_phi < 360; cell_phi++){
            for(int cell_theta = min(lower_theta,360-upper_theta); cell_theta < 180; cell_theta++){
                if(grid[cell_phi*180+cell_theta].cnt<=0)
                {
                    continue;
                }

                if(calc_angle(current_satellite,grid[cell_phi*180+cell_theta].cells[0])>current_satellite->delta*1.5)
                {
                    continue;
                }

                proc_sector(grid+cell_phi*180+cell_theta,current_satellite,satellites_end,satellitesTotalCount);
                }
        }

        return 1;
    }
    //Если спутник покрывает южный полюс
    if(lower_theta < 0){
        for(int cell_phi = 0; cell_phi < 360; cell_phi++){
            for(int cell_theta = 0; cell_theta < max(-lower_theta,upper_theta); cell_theta++){
                if(grid[cell_phi*180+cell_theta].cnt<=0)
                {
                    continue;
                }

                if(calc_angle(current_satellite,grid[cell_phi*180+cell_theta].cells[0])>current_satellite->delta*1.5)
                {
                    continue;
                }

                proc_sector(grid+cell_phi*180+cell_theta,current_satellite,satellites_end,satellitesTotalCount);
                }
        }

        return 1;
    }
    //Обычный случай
    for(int i = lower_phi; i < upper_phi; i++)
    {
        int cell_phi = i%360;
        if(cell_phi < 0)
        {
            cell_phi+=360;
        }
        for(int cell_theta = lower_theta; cell_theta < upper_theta; cell_theta++)
        {
            if(grid[cell_phi*180+cell_theta].cnt<=0)
            {
                continue;
            }

            if(calc_angle(current_satellite,grid[cell_phi*180+cell_theta].cells[0])>current_satellite->delta*1.5)
            {
                continue;
            }

            proc_sector(grid+cell_phi*180+cell_theta,current_satellite,satellites_end,satellitesTotalCount);
        }
    }
    return 1;
}

int main()
{


    auto start = chrono::high_resolution_clock::now();

    //Инициализация карты
    grid = (HSECTOR)calloc(360*180,sizeof(SECTOR));

    //Считывание клеток
    ifstream cells_file("cells.txt");
    //cells_file.tie(0);
    cells_file >> cellsTotalCount;
    cells = (HCELL)malloc(cellsTotalCount*sizeof(CELL));
    for (HCELL cell = cells; cell < cells + cellsTotalCount; cell++)
    {
        cells_file >> cell->x >> cell->y >> cell->z;
    }
    cells_file.close();

    for (HCELL cell = cells; cell < cells + cellsTotalCount; cell++)
    {
        proc_Cell(cell,grid);
    }
    //Считывание спутников
    ifstream satellites_file("satellites.txt");
    //satellites_file.tie(0);
    satellites_file >> satellitesTotalCount >> timesTotalCount;
    satellites = (HSATELLITE)malloc(satellitesTotalCount*timesTotalCount*sizeof(SATELLITE));
    for(HSATELLITE satellite = satellites; satellite < satellites + satellitesTotalCount*timesTotalCount; satellite++)
    {
        satellites_file >> satellite->x >> satellite->y >> satellite->z;
        proc_Satellite(satellite);
    }
    //cout << satellites[872].phi <<" " << satellites[872].theta <<" " << satellites[872].delta <<endl;
    satellites_file.close();
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
    cout << "Reading files finished in "<< duration.count() / 1000.0 / 1000.0 / 1000.0<< " seconds." << endl;
    start = chrono::high_resolution_clock::now();
    ofstream result("result.txt");
    for(int t = 0; t < timesTotalCount; t++)
    {
        if(t%100 == 0){
            cout << "Epoch:" << t<< endl;
        }

        for(HSATELLITE satellite = satellites+t*satellitesTotalCount; satellite < satellites+(t+1)*satellitesTotalCount; satellite++)
        {
            proc_time(satellite,satellites+satellitesTotalCount*timesTotalCount,satellitesTotalCount,grid);
        }
        //Процедура заполнения файла номерами спутников, потдерживающих данную клетку
        for(HCELL cell = cells; cell < cells + cellsTotalCount; cell++){
            if(cell->sattelite_id == -1 && cell->max_satellite >=0){
                cell->sattelite_id = cell->max_satellite;
            }
            result << cell->sattelite_id << " ";
            if(cell->sattelite_id>=0){
                cell->connection_time--;
            }
            if(cell->connection_time<=0){
                cell->connection_time = 0;
                cell->max_satellite = -1;
                cell->sattelite_id = -1;
            }
        }
        result << endl;
    }
    result.close();
    end = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::nanoseconds>(end - start);
    cout << "Processing satellites and writing finished in: "<< duration.count() / 1000.0 / 1000.0 / 1000.0<< " seconds." << endl;

    return 0;
}
