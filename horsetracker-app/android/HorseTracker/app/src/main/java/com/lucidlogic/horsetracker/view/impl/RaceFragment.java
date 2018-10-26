package com.lucidlogic.horsetracker.view.impl;

import android.content.Context;
import android.net.Uri;
import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.support.v7.widget.LinearLayoutManager;
import android.support.v7.widget.RecyclerView;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;

import com.lucidlogic.horsetracker.R;
import com.lucidlogic.horsetracker.adapter.RunnerRecyclerViewAdapter;
import com.lucidlogic.horsetracker.model.Entity;
import com.lucidlogic.horsetracker.model.Race;
import com.lucidlogic.horsetracker.model.Runner;
import com.lucidlogic.horsetracker.presenter.RacePresenter;
import com.lucidlogic.horsetracker.presenter.RacecardPresenter;
import com.lucidlogic.horsetracker.presenter.impl.RacePresenterImpl;
import com.lucidlogic.horsetracker.presenter.impl.RacecardPresenterImpl;
import com.lucidlogic.horsetracker.view.RaceView;

public class RaceFragment extends Fragment implements RaceView {

    private static final String RACE_OBJECT_ID = "RACE_OBJECT_ID";

    private RecyclerView mRecyclerView;
    private RecyclerView.Adapter mAdapter;
    private RecyclerView.LayoutManager mLayoutManager;

    private RacePresenter racePresenter;

    private String raceObjectId;

    public RaceFragment() {
        racePresenter = new RacePresenterImpl();
        racePresenter.setView(this);
    }

    public static RaceFragment newInstance(String raceObjectId) {
        RaceFragment fragment = new RaceFragment();
        Bundle args = new Bundle();
        args.putString(RACE_OBJECT_ID, raceObjectId);
        fragment.setArguments(args);
        return fragment;
    }

    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        if (getArguments() != null) {
            raceObjectId = getArguments().getString(RACE_OBJECT_ID);
        }
    }

    @Override
    public View onCreateView(LayoutInflater inflater, ViewGroup container,
                             Bundle savedInstanceState) {
        View v = inflater.inflate(R.layout.fragment_race_list, container, false);
        mRecyclerView = (RecyclerView) v.findViewById(R.id.runner_list);
        racePresenter.updateView(raceObjectId);
        return v;
    }


    @Override
    public void onDetach() {
        super.onDetach();
        racePresenter = null;
    }

    @Override
    public void onDestroy() {
        super.onDestroy();
    }

    @Override
    public void updateView(Race race) {
        mLayoutManager = new LinearLayoutManager(this.getActivity());
        mRecyclerView.setLayoutManager(mLayoutManager);
        mAdapter = new RunnerRecyclerViewAdapter(race.getRunners());
        mRecyclerView.setAdapter(mAdapter);
    }
}
