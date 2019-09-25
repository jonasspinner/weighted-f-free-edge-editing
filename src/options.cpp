//
// Created by jonas on 24.09.19.
//

#include "options.h"


namespace Options {

    std::istream& operator>>(std::istream& in, Selector& selector) {
        std::string token;
        in >> token;
        if (token == "FirstFound")
            selector = Selector::FirstFound;
        else if (token == "LeastWeight")
            selector = Selector::LeastWeight;
        else if (token == "MostMarkedPairs")
            selector = Selector::MostMarkedPairs;
        else
            in.setstate(std::ios_base::failbit);
        return in;
    }

    std::ostream &operator<<(std::ostream &os, Selector selector) {
        switch (selector) {
            case Selector::FirstFound:
                return os << "FirstFound";
            case Selector::LeastWeight:
                return os << "LeastWeight";
            case Selector::MostMarkedPairs:
                return os << "MostMarkedPairs";
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
        else if (token == "C4P4")
            fsg = FSG::C4P4;
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
            case FSG::C4P4:
                return os << "C4P4";
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
        if (token == "Trivial")
            lower_bound = LB::Trivial;
        else if (token == "Greedy")
            lower_bound = LB::Greedy;
        else if (token == "SortedGreedy")
            lower_bound = LB::SortedGreedy;
        else if (token == "LocalSearch")
            lower_bound = LB::LocalSearch;
        else if (token == "LPRelaxation")
            lower_bound = LB::LPRelaxation;
        else
            in.setstate(std::ios_base::failbit);
        return in;
    }

    std::ostream &operator<<(std::ostream &os, LB lower_bound) {
        switch (lower_bound) {
            case LB::Trivial:
                return os << "Trivial";
            case LB::Greedy:
                return os << "Greedy";
            case LB::SortedGreedy:
                return os << "SortedGreedy";
            case LB::LocalSearch:
                return os << "LocalSearch";
            case LB::LPRelaxation:
                return os << "LPRelaxation";
            default:
                return os;
        }
    }

    YAML::Emitter &operator<<(YAML::Emitter &out, LB lower_bound) {
        std::ostringstream ss;
        ss << lower_bound;
        return out << ss.str();
    }

    std::istream &operator>>(std::istream &in, SolverType &type) {
        std::string token;
        in >> token;
        if (token == "FPT")
            type = SolverType::FPT;
        else if (token == "ILP")
            type = SolverType::ILP;
        else
            in.setstate(std::ios_base::failbit);
        return in;
    }

    std::ostream &operator<<(std::ostream &os, SolverType type) {
        switch (type) {
            case SolverType::FPT:
                return os << "FPT";
            case SolverType::ILP:
                return os << "ILP";
            default:
                return os;
        }
    }

    YAML::Emitter &operator<<(YAML::Emitter &out, SolverType type) {
        std::ostringstream ss;
        ss << type;
        return out << ss.str();
    }
}
