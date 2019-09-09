//
// Created by jonas on 17.07.19.
//


#include "Configuration.h"


namespace Options {

    std::istream& operator>>(std::istream& in, Selector& selector) {
        std::string token;
        in >> token;
        if (token == "FirstEditable")
            selector = Selector::FirstEditable;
        else if (token == "LeastWeight")
            selector = Selector::LeastWeight;
        else
            in.setstate(std::ios_base::failbit);
        return in;
    }

    std::ostream &operator<<(std::ostream &os, Selector selector) {
        switch (selector) {
            case Selector::FirstEditable:
                return os << "FirstEditable";
            case Selector::LeastWeight:
                return os << "LeastWeight";
            default:
                return os;
        }
    }

    YAML::Emitter &operator<<(YAML::Emitter &out, Selector selector) {
        std::ostringstream ss;
        ss << selector;
        return out << ss.str();
    }
}


namespace Options {

    std::istream& operator>>(std::istream& in, FSG& fsg) {
        std::string token;
        in >> token;
        if (token == "P3")
            fsg = FSG::P3;
        else if (token == "P4C4")
            fsg = FSG::P4C4;
        else if (token == "C4_C5_2K2")
            fsg = FSG::C4_C5_2K2;
        else if (token == "C4_C5_P5_Bowtie_Necktie")
            fsg = FSG::C4_C5_P5_Bowtie_Necktie;
        else
            in.setstate(std::ios_base::failbit);
        return in;
    }

    std::ostream &operator<<(std::ostream &os, FSG fsg) {
        switch (fsg) {
            case FSG::P3:
                return os << "P3";
            case FSG::P4C4:
                return os << "P4C4";
            case FSG::C4_C5_2K2:
                return os << "C4_C5_2K2";
            case FSG::C4_C5_P5_Bowtie_Necktie:
                return os << "C4_C5_P5_Bowtie_Necktie";
            default:
                return os;
        }
    }

    YAML::Emitter &operator<<(YAML::Emitter &out, FSG fsg) {
        std::ostringstream ss;
        ss << fsg;
        return out << ss.str();
    }
}


namespace Options {

    std::istream& operator>>(std::istream& in, LB& lower_bound) {
        std::string token;
        in >> token;
        if (token == "No")
            lower_bound = LB::No;
        else if (token == "Greedy")
            lower_bound = LB::Greedy;
        else if (token == "LocalSearch")
            lower_bound = LB::LocalSearch;
        else if (token == "LinearProgram")
            lower_bound = LB::LinearProgram;
        else
            in.setstate(std::ios_base::failbit);
        return in;
    }

    std::ostream &operator<<(std::ostream &os, LB lower_bound) {
        switch (lower_bound) {
            case LB::No:
                return os << "No";
            case LB::Greedy:
                return os << "Greedy";
            case LB::LocalSearch:
                return os << "LocalSearch";
            case LB::LinearProgram:
                return os << "LinearProgram";
            default:
                return os;
        }
    }

    YAML::Emitter &operator<<(YAML::Emitter &out, LB lower_bound) {
        std::ostringstream ss;
        ss << lower_bound;
        return out << ss.str();
    }
}
