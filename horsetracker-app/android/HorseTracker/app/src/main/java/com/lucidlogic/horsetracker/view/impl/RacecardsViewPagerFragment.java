package com.lucidlogic.horsetracker.view.impl;


import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;

import com.lucidlogic.horsetracker.R;

public class RacecardsViewPagerFragment extends Fragment {

    public RacecardsViewPagerFragment() {
    }


    public static RacecardsViewPagerFragment newInstance() {
        return new RacecardsViewPagerFragment();
    }

    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);

    }

    @Override
    public View onCreateView(LayoutInflater inflater, ViewGroup container,
                             Bundle savedInstanceState) {
        return inflater.inflate(R.layout.fragment_racecards_view_pager, container, false);
    }

}
