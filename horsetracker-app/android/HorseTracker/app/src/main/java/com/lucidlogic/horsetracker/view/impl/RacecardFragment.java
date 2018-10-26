package com.lucidlogic.horsetracker.view.impl;

import android.content.Context;
import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.widget.ExpandableListView;

import com.lucidlogic.horsetracker.R;
import com.lucidlogic.horsetracker.adapter.RacecardExpandableListAdapter;
import com.lucidlogic.horsetracker.model.Racecard;
import com.lucidlogic.horsetracker.presenter.RacecardPresenter;
import com.lucidlogic.horsetracker.presenter.impl.RacecardPresenterImpl;
import com.lucidlogic.horsetracker.view.RacecardView;

public class RacecardFragment extends Fragment implements RacecardView {

    RacecardPresenter racecardPresenter;
    ExpandableListView expListView;

    public RacecardFragment() {
        racecardPresenter = new RacecardPresenterImpl();
        racecardPresenter.setView(this);
    }

    public static RacecardFragment newInstance() {
        return new RacecardFragment();
    }

    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
    }

    @Override
    public View onCreateView(LayoutInflater inflater, ViewGroup container,
                             Bundle savedInstanceState) {
        View view = inflater.inflate(R.layout.fragment_racecard_list, container, false);
        expListView = (ExpandableListView) view.findViewById(R.id.lvExp);
        racecardPresenter.updateView();
        return view;
    }


    @Override
    public void onDetach() {
        super.onDetach();
        racecardPresenter.setView(null);
    }

    @Override
    public void updateView(Racecard racecard) {
        expListView.setAdapter(new RacecardExpandableListAdapter(this.getContext(), racecard));
    }

}
