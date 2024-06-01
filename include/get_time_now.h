#ifndef GET_TIME_NOW_H
#define GET_TIME_NOW_H

#include <string>
#include <chrono>

inline std::string get_time_now()
{
    auto now = std::chrono::system_clock::now();
    time_t tt = std::chrono::system_clock::to_time_t(now);
    struct tm *localTime = localtime(&tt);
    char time_buffer_str[80];

    strftime(time_buffer_str, sizeof(time_buffer_str), "%Y-%m-%d %H:%M:%S", localTime);

    return std::string(time_buffer_str);
}

#endif // GET_TIME_NOW_H